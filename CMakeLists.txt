#
# Copyright 2021-2023 Software Radio Systems Limited
#
# This file is part of srsRAN
#
# srsRAN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# srsRAN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# A copy of the GNU Affero General Public License can be found in
# the LICENSE file in the top-level directory of this distribution
# and at http://www.gnu.org/licenses/.
#


add_executable(ldpc_encoder_benchmark ldpc_encoder_benchmark.cpp)

target_include_directories(ldpc_encoder_benchmark PRIVATE /home/ubuntu/linux/tools/lib/perf/include)
target_link_libraries(ldpc_encoder_benchmark srsran_channel_coding srslog /home/ubuntu/perf_cycle_codes/pmuv3_plugin.a /home/ubuntu/perf_cycle_codes/libperf.a /home/ubuntu/perf_cycle_codes/libapi.a)


add_test(ldpc_encoder_benchmark ldpc_encoder_benchmark -s -R 1 -T generic)

add_executable(ldpc_decoder_benchmark ldpc_decoder_benchmark.cpp)
target_link_libraries(ldpc_decoder_benchmark srsran_channel_coding srslog)
add_test(ldpc_decoder_benchmark ldpc_decoder_benchmark -s -R 1 -I 2 -T generic)


