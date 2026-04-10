/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace ijedi {

extern "C" {

void ijedi_mpas_state_read_field_f90(const void * geom,
                                     const char * filepath,
                                     const char * stream_name,
                                     const char * var_name,
                                     const int & nvals,
                                     const int * global_indices,
                                     const int & nlevels,
                                     double * values,
                                     const double & scaling);

void ijedi_mpas_state_write_field_f90(const void * geom,
                                      const char * filepath,
                                      const char * stream_name,
                                      const char * var_name,
                                      const int & nvals,
                                      const int * global_indices,
                                      const int & nlevels,
                                      const double * values);

void ijedi_mpas_state_write_flush_f90(const void * geom,
                                      const char * filepath,
                                      const char * stream_name);

}  // extern "C"

}  // namespace ijedi
