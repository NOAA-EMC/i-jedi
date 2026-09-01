"""NeuralGCM (Google Research) as an i-jedi forecast model.

Everything model-specific lives here, which is the point of the adapter layer: i-jedi's C++
knows nothing about NeuralGCM, and adding AIGFS or a PyTorch model means writing a sibling of
this file, not touching C++.

Array conventions, which is where the real work is
--------------------------------------------------
i-jedi hands over flat float64 buffers, level-major over atlas's global node ordering. For the
Gaussian grids these models use, atlas numbers latitudes from north to south, so

    buffer.reshape(n_levels, n_lat, n_lon)

has latitude descending. NeuralGCM wants latitude ascending and longitude before latitude, so
each field is flipped along latitude and transposed. Getting either wrong produces a field
that still looks plausible - a hemisphere-swapped temperature field is smooth and has the
right range - which is why the test compares against a rollout that never went through this
code rather than merely checking the forecast completes.

The grid itself needs no interpolation: NeuralGCM's TL63 grid and atlas's
"regular_gaussian N=32" agree to 3e-14 degrees.
"""

from __future__ import annotations

import datetime
from typing import Dict

import numpy as np

from ijedi_model_adapter import Declaration, IjediModelAdapter, register

# NeuralGCM's variable names -> the JEDI long names they correspond to. i-jedi uses these as
# defaults and the yaml can override any of them.
INPUT_VARIABLES = {
    "u_component_of_wind": "eastward_wind",
    "v_component_of_wind": "northward_wind",
    "temperature": "air_temperature",
    "geopotential": "geopotential",
    "specific_humidity": "water_vapor_mixing_ratio_wrt_moist_air",
    "specific_cloud_liquid_water_content": "cloud_liquid_water",
    "specific_cloud_ice_water_content": "cloud_liquid_ice",
}
FORCING_VARIABLES = {
    "sea_surface_temperature": "sea_surface_temperature",
    "sea_ice_cover": "sea_ice_area_fraction",
}


@register("neuralgcm")
class NeuralGcmAdapter(IjediModelAdapter):
    """Configuration keys, all under the yaml's "adapter config":

        checkpoint   path to a pickled checkpoint, or "demo" for the TL63 model packaged
                     with the neuralgcm wheel (needs no network, used by the tests)
        rng seed     integer seed for the stochastic model. Fixed and logged, because a
                     forecast whose noise draw is not recorded is not reproducible, and
                     because the adjoint must later linearise about the same draw.
    """

    def __init__(self, config):
        super().__init__(config)

        # Imported here, not at module scope, so that listing the available adapters does
        # not require JAX to be installed.
        import jax
        import neuralgcm

        self._jax = jax
        self._np = np

        checkpoint = config.get("checkpoint", "demo")
        if checkpoint == "demo":
            from neuralgcm import demo
            ckpt = demo.load_checkpoint_tl63_stochastic()
            self._source = "packaged demo TL63 checkpoint"
        else:
            import pickle
            with open(checkpoint, "rb") as handle:
                ckpt = pickle.load(handle)
            self._source = checkpoint

        self.model = neuralgcm.PressureLevelModel.from_checkpoint(ckpt)
        self.seed = int(config.get("rng seed", 0))

        coords = self.model.data_coords
        self.lat_deg = np.degrees(np.asarray(coords.horizontal.latitudes))
        self.lon_deg = np.degrees(np.asarray(coords.horizontal.longitudes))
        self.levels_hpa = np.asarray(coords.vertical.centers, dtype=float)
        self.n_lat = len(self.lat_deg)
        self.n_lon = len(self.lon_deg)
        self.n_lev = len(self.levels_hpa)

        # NeuralGCM reports its timestep as a numpy timedelta64.
        self.timestep_seconds = int(
            np.asarray(self.model.timestep).astype("timedelta64[s]").astype(int))

        self._state = None
        self._forcings = None
        self._forcings_series = None
        self._trajectory = None
        self._trajectory_time = None
        self._rng = None
        self._sim_time = None

        import logging
        logging.getLogger(__name__).info(
            "NeuralGCM adapter: %s, %d levels, %dx%d grid, timestep %ds, seed %d",
            self._source, self.n_lev, self.n_lon, self.n_lat,
            self.timestep_seconds, self.seed)

    # -- declaration ------------------------------------------------------------------------

    @property
    def declaration(self) -> Declaration:
        return Declaration(
            description=f"NeuralGCM ({self._source})",
            timestep=datetime.timedelta(seconds=self.timestep_seconds),
            # i-jedi works in Pa; NeuralGCM labels its levels in hPa.
            levels_pa=[float(v) * 100.0 for v in self.levels_hpa],
            input_variables=dict(INPUT_VARIABLES),
            forcing_variables=dict(FORCING_VARIABLES),
            supports_linear=True,
        )

    # -- layout conversion ------------------------------------------------------------------

    def _to_model(self, flat: np.ndarray, levels: int) -> np.ndarray:
        """i-jedi's (level, lat descending, lon) -> NeuralGCM's (level, lon, lat ascending)."""
        if levels > 1:
            a = flat.reshape(levels, self.n_lat, self.n_lon)
            return a[:, ::-1, :].transpose(0, 2, 1)
        a = flat.reshape(self.n_lat, self.n_lon)
        return a[::-1, :].transpose(1, 0)

    def _from_model(self, arr: np.ndarray) -> np.ndarray:
        """NeuralGCM's layout back to i-jedi's, flattened."""
        a = np.asarray(arr)
        if a.ndim == 3:
            return np.ascontiguousarray(a.transpose(0, 2, 1)[:, ::-1, :]).ravel()
        return np.ascontiguousarray(a.transpose(1, 0)[::-1, :]).ravel()

    def _to_model_jax(self, flat, levels: int):
        """As _to_model, but written so JAX can trace through it."""
        import jax.numpy as jnp
        if levels > 1:
            a = jnp.reshape(flat, (levels, self.n_lat, self.n_lon))
            return jnp.transpose(a[:, ::-1, :], (0, 2, 1))
        a = jnp.reshape(flat, (self.n_lat, self.n_lon))
        return jnp.transpose(a[::-1, :], (1, 0))

    def _from_model_jax(self, arr):
        """As _from_model, but traceable."""
        import jax.numpy as jnp
        a = jnp.asarray(arr)
        if a.ndim == 3:
            return jnp.reshape(jnp.transpose(a, (0, 2, 1))[:, ::-1, :], (-1,))
        return jnp.reshape(jnp.transpose(a, (1, 0))[::-1, :], (-1,))

    def _dataset(self, fields: Dict[str, np.ndarray], valid_time: datetime.datetime):
        """Build the xarray NeuralGCM's own converters expect.

        A time dimension is required even though only one time is present: the forcings are
        interpolated in time during a rollout, and a dataset without a time axis produces
        0-d forcing arrays and an obscure failure inside jnp.interp.
        """
        import xarray

        data = {}
        for name, flat in fields.items():
            if name in INPUT_VARIABLES:
                data[name] = (("time", "level", "longitude", "latitude"),
                              self._to_model(flat, self.n_lev)[np.newaxis])
            elif name in FORCING_VARIABLES:
                data[name] = (("time", "longitude", "latitude"),
                              self._to_model(flat, 1)[np.newaxis])
        return xarray.Dataset(
            data,
            coords={"time": np.array([np.datetime64(valid_time.replace(tzinfo=None), "ns")]),
                    "level": self.levels_hpa,
                    "longitude": self.lon_deg,
                    "latitude": self.lat_deg},
        )

    # -- nonlinear --------------------------------------------------------------------------

    def encode(self, inputs, valid_time: datetime.datetime) -> None:
        ds = self._dataset(inputs, valid_time)
        self._trajectory_time = np.array(
            [np.datetime64(valid_time.replace(tzinfo=None), "ns")])
        model_inputs = self.model.inputs_from_xarray(ds.isel(time=0))
        # encode and decode take the forcings at one time; unroll takes the series, since it
        # interpolates them to each internal step.
        self._forcings = self.model.forcings_from_xarray(ds.isel(time=0))
        self._forcings_series = self.model.forcings_from_xarray(ds)
        self._rng = self._jax.random.key(self.seed)
        # The model carries its own notion of simulation time alongside the fields; the
        # traced linear map reuses the value the trajectory was built with.
        self._sim_time = model_inputs["sim_time"]
        self._state = self.model.encode(model_inputs, self._forcings, self._rng)

    def advance(self, steps: int) -> None:
        if self._state is None:
            raise RuntimeError("advance called before encode")
        if steps <= 0:
            return
        # One unroll of the whole interval rather than a Python loop over single steps: the
        # loop stays inside the compiled region, and asking for a single output avoids
        # materialising the intermediate states.
        timedelta = np.timedelta64(self.timestep_seconds * int(steps), "s")
        self._state, _ = self.model.unroll(
            self._state, self._forcings_series, steps=1, timedelta=timedelta,
            start_with_input=False)

    def decode(self, outputs) -> None:
        if self._state is None:
            raise RuntimeError("decode called before encode")
        # unroll returns a leading singleton output axis; drop it.
        state = self._jax.tree_util.tree_map(
            lambda x: x[0] if getattr(x, "ndim", 0) and x.shape[0] == 1 else x, self._state)
        decoded = self.model.decode(state, self._forcings)
        ds = self.model.data_to_xarray(decoded, times=None)
        for name in INPUT_VARIABLES:
            if name in ds:
                outputs[name] = self._from_model(ds[name].values)

    def reset(self) -> None:
        self._state = None
        self._forcings = None
        self._forcings_series = None
        self._trajectory = None

    # -- tangent linear and adjoint ---------------------------------------------------------
    #
    # The operator being linearised is the whole composition decode(advance^n(encode(x))),
    # not advance alone. encode and decode are learned networks in NeuralGCM, so leaving them
    # out would linearise something the assimilation never actually applies.
    #
    # The random draw is pinned to the same seed the trajectory used. For a stochastic
    # checkpoint this matters and is easy to get wrong silently: jvp and vjp are exact
    # transposes of whatever function they are handed, so an adjoint taken about a different
    # noise realisation still passes a dot-product test while being the adjoint of the wrong
    # operator.

    def _forward(self, fields_flat, steps: int):
        """The full pressure-level to pressure-level map, as a function JAX can differentiate.

        Takes and returns a dict of flat i-jedi-ordered arrays, so the perturbation the
        assimilation carries lives in the same space as the state it perturbs.

        Note this deliberately avoids inputs_from_xarray and data_to_xarray, which the
        nonlinear path uses happily. Those route through xarray, and reading a DataArray's
        .values calls np.asarray, which raises TracerArrayConversionError on a JAX tracer.
        The model's underlying interface is a plain dict of (level, longitude, latitude)
        arrays plus a sim_time scalar, so the traced path builds that directly.
        """
        inputs = {name: self._to_model_jax(fields_flat[name], self.n_lev).astype(
                      self._jax.numpy.float32)
                  for name in INPUT_VARIABLES if name in fields_flat}
        inputs["sim_time"] = self._sim_time

        state = self.model.encode(inputs, self._forcings, self._rng)
        if steps > 0:
            timedelta = np.timedelta64(self.timestep_seconds * int(steps), "s")
            state, _ = self.model.unroll(
                state, self._forcings_series, steps=1, timedelta=timedelta,
                start_with_input=False)
            state = self._jax.tree_util.tree_map(
                lambda x: x[0] if getattr(x, "ndim", 0) and x.shape[0] == 1 else x, state)
        decoded = self.model.decode(state, self._forcings)
        return {name: self._from_model_jax(decoded[name])
                for name in INPUT_VARIABLES if name in decoded}

    def set_trajectory(self, inputs, valid_time: datetime.datetime) -> None:
        """Remember the point to linearise about, and the forcings and noise that go with it."""
        self.encode(inputs, valid_time)
        self._rng = self._jax.random.key(self.seed)
        self._trajectory = {k: np.asarray(v, dtype=np.float64).copy()
                            for k, v in inputs.items() if k in INPUT_VARIABLES}

    def advance_tl(self, dx_in, dx_out, steps: int) -> None:
        if self._trajectory is None:
            raise RuntimeError("advance_tl called before set_trajectory")
        primals = {k: self._jax.numpy.asarray(v) for k, v in self._trajectory.items()}
        tangents = {k: self._jax.numpy.asarray(dx_in[k]) for k in primals}
        _, out = self._jax.jvp(lambda x: self._forward(x, steps), (primals,), (tangents,))
        for name, value in out.items():
            dx_out[name] = np.asarray(value, dtype=np.float64).ravel()

    def advance_ad(self, dx_out, dx_in, steps: int) -> None:
        if self._trajectory is None:
            raise RuntimeError("advance_ad called before set_trajectory")
        primals = {k: self._jax.numpy.asarray(v) for k, v in self._trajectory.items()}
        _, vjp = self._jax.vjp(lambda x: self._forward(x, steps), primals)
        seed = {k: self._jax.numpy.asarray(dx_out[k]) for k in primals}
        (grad,) = vjp(seed)
        for name, value in grad.items():
            dx_in[name] = np.asarray(value, dtype=np.float64).ravel()
