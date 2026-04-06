/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace eckit
{
    namespace mpi
    {
        class Comm;
    }  // namespace mpi
    class Configuration;
}  // namespace eckit

namespace ijedi
{
    extern "C"
    {
        void ijedi_mpas_geom_setup_f90(void *& geom, const eckit::Configuration &,
			               const eckit::mpi::Comm *);
        void ijedi_mpas_geom_clone_f90(void *& geom, const void * other);
        void ijedi_mpas_geom_delete_f90(void *& geom);
        void ijedi_mpas_geom_get_num_nodes_and_elements_f90(const void * geom, int & num_nodes,
                                       int & num_tris);
        void ijedi_mpas_geom_get_coords_and_connectivities_f90(const void * geom, 
			               const int & num_nodes,
                                       double * lons, double * lats, int * ghosts,
                                       int * global_indices, int * remote_indices,
                                       int * partition, const int & num_tri_nodes,
                                       int * raw_tri_boundary_nodes);
        void ijedi_mpas_geom_get_vertical_resolution_f90(const void * geom, int & nVertLevels);
        void ijedi_mpas_geom_get_local_cell_counts_f90(const void * geom, int & nCells, i
			               int & nCellsSolve);
        void ijedi_mpas_geom_get_global_cell_count_f90(const void * geom, int & nCellsGlobal);
        void ijedi_mpas_geom_get_area_f90(const void * geom, const int & n, double * area);
    }  // extern "C"

}  // namespace ijedi
