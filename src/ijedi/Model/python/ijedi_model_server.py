#!/usr/bin/env python3
"""Out-of-process model server: drives an adapter on behalf of an i-jedi process.

i-jedi spawns this with one end of a socketpair on file descriptor 3 and talks to it in a
simple request/response protocol. A socket rather than the standard streams, because JAX,
absl and friends write to stdout and stderr whenever they feel like it and would corrupt a
protocol that shared them. Here stdout and stderr stay free for exactly that chatter, and are
inherited from the parent so it shows up in the JEDI log where someone will see it.

Wire format, both directions:

    8 bytes   little-endian unsigned length of the JSON header
    N bytes   JSON header, UTF-8
    M bytes   payload: float64 arrays, little-endian, concatenated in the order the
              header's "arrays" list gives, each of the stated size

Requests carry {"op": ..., ...}; responses carry {"status": "ok"} plus whatever the operation
returns, or {"status": "error", "message": ..., "traceback": ...}. An error is a normal
response, not a dropped connection: the C++ side turns it into an exception with the Python
traceback attached, which is far more use than "the model server died".
"""

from __future__ import annotations

import json
import os
import socket
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ijedi_model_host import Host, error_response  # noqa: E402

_LEN = struct.Struct("<Q")
CONTROL_FD = 3


def _recv_exactly(sock: socket.socket, n: int) -> bytes:
    chunks = []
    got = 0
    while got < n:
        chunk = sock.recv(min(1 << 20, n - got))
        if not chunk:
            raise ConnectionError("i-jedi closed the connection")
        chunks.append(chunk)
        got += len(chunk)
    return b"".join(chunks)


def read_message(sock: socket.socket):
    (length,) = _LEN.unpack(_recv_exactly(sock, _LEN.size))
    header = json.loads(_recv_exactly(sock, length).decode("utf-8"))

    arrays = {}
    for spec in header.get("arrays", []):
        buf = _recv_exactly(sock, int(spec["size"]) * 8)
        arrays[spec["name"]] = np.frombuffer(buf, dtype="<f8").copy()
    return header, arrays


def write_message(sock: socket.socket, header: dict, arrays=None) -> None:
    arrays = arrays or {}
    if arrays:
        header = dict(header)
        header["arrays"] = [{"name": k, "size": int(v.size)} for k, v in arrays.items()]
    blob = json.dumps(header).encode("utf-8")
    sock.sendall(_LEN.pack(len(blob)) + blob)
    for spec in header.get("arrays", []):
        sock.sendall(np.ascontiguousarray(arrays[spec["name"]], dtype="<f8").tobytes())


def main() -> int:
    try:
        sock = socket.socket(fileno=os.dup(CONTROL_FD))
    except OSError as exc:
        print(f"ijedi_model_server: no control socket on fd {CONTROL_FD} ({exc})",
              file=sys.stderr)
        return 1

    host = None
    while True:
        try:
            header, arrays = read_message(sock)
        except (ConnectionError, OSError):
            return 0  # parent went away; nothing to report

        op = header.get("op")
        if op == "shutdown":
            try:
                write_message(sock, {"status": "ok"})
            except OSError:
                pass
            return 0

        try:
            if op == "init":
                host = Host(header["adapter"], header.get("config", {}))
                write_message(sock, {"status": "ok"})
                continue
            if host is None:
                raise RuntimeError(f"received '{op}' before 'init'")

            response = host.call(op, header, arrays)
            out_arrays = response.pop("arrays", None)
            write_message(sock, response, out_arrays)
        except BaseException as exc:  # noqa: BLE001 - every failure must reach the caller
            try:
                write_message(sock, error_response(exc))
            except OSError:
                return 1


if __name__ == "__main__":
    sys.exit(main())
