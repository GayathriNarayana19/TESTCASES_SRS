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
#link_directories(/home/ubuntu/perf_cycle_codes/)
link_directories(/home/ubuntu/linux-5.17.4/tools/lib/perf/ /home/ubuntu/linux-5.17.4/tools/lib/api/ /home/ubuntu/perf_cycle_codes/)
target_link_libraries(ldpc_encoder_benchmark srsran_channel_coding srslog /home/ubuntu/perf_cycle_codes/pmuv3_dummy.a /home/ubuntu/linux-5.17.4/tools/lib/perf/libperf.a /home/ubuntu/linux-5.17.4/tools/lib/api/libapi.a)
#target_link_libraries(ldpc_encoder_benchmark srsran_channel_coding srslog /home/ubuntu/perf_cycle_codes/pmuv3.a /home/ubuntu/linux-5.17.4/tools/lib/perf/libperf.a /home/ubuntu/linux-5.17.4/tools/lib/api/libapi.a)


target_link_libraries(ldpc_encoder_benchmark srsran_channel_coding srslog /home/ubuntu/perf_cycle_codes/pmuv3_plugin.a /home/ubuntu/linux-5.17.4/tools/lib/perf/libperf.a /home/ubuntu/linux-5.17.4/tools/lib/api/libapi.a)

#target_link_libraries(ldpc_encoder_benchmark srsran_channel_coding srslog /home/ubuntu/linux-5.17.4/tools/lib/perf/libperf.a /home/ubuntu/linux-5.17.4/tools/lib/api/libapi.a /home/ubuntu/perf_cycle_codes/pmuv3.a)
#target_link_libraries(ldpc_encoder_benchmark /home/ubuntu/linux-5.17.4/tools/lib/perf/libperf.a) #added by Gayathri
#target_link_libraries(ldpc_encoder_benchmark /home/ubuntu/linux-5.17.4/tools/lib/api/libapi.a) #added by Gayathri
target_include_directories(ldpc_encoder_benchmark PRIVATE /home/ubuntu/linux-5.17.4/tools/)
target_include_directories(ldpc_encoder_benchmark PRIVATE /home/ubuntu/linux-5.17.4/tools/lib/perf/include/perf)
#include_directories(/home/ubuntu/linux-5.17.4/include)
add_test(ldpc_encoder_benchmark ldpc_encoder_benchmark -s -R 1 -T generic)

add_executable(ldpc_decoder_benchmark ldpc_decoder_benchmark.cpp)
target_link_libraries(ldpc_decoder_benchmark srsran_channel_coding srslog)
add_test(ldpc_decoder_benchmark ldpc_decoder_benchmark -s -R 1 -I 2 -T generic)

add_executable(testcase testcase.c)
#include_directories(/home/ubuntu/linux-5.17.4/tools/lib/perf)
target_link_libraries(testcase /home/ubuntu/linux-5.17.4/tools/lib/perf/libperf.a)
target_link_libraries(testcase /home/ubuntu/linux-5.17.4/tools/lib/api/libapi.a)
add_executable(pmuv3_cycle pmuv3_cycle.c)
#include_directories(/home/ubuntu/linux-5.17.4/tools/lib/perf)
target_link_libraries(pmuv3_cycle /home/ubuntu/linux-5.17.4/tools/lib/perf/libperf.a)
target_link_libraries(pmuv3_cycle /home/ubuntu/linux-5.17.4/tools/lib/api/libapi.a)
set(CMAKE_VERBOSE_MAKEFILE ON)

