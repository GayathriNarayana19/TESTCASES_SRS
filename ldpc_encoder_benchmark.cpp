/*
 *
 * Copyright 2021-2023 Software Radio Systems Limited
 *
 * This file is part of srsRAN.
 *
 * srsRAN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * srsRAN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * A copy of the GNU Affero General Public License can be found in
 * the LICENSE file in the top-level directory of this distribution
 * and at http://www.gnu.org/licenses/.
 *
 */

#include "srsran/phy/upper/channel_coding/channel_coding_factories.h"
#include "srsran/support/benchmark_utils.h"
#include "srsran/support/srsran_test.h"
#include <getopt.h>
#include <random>
#include <iostream>
#include <cstdint>
#include <climits>
#include <chrono>
#include <vector>
#include <algorithm>
#include <fstream> 
#include <iomanip>
//#include </home/ubuntu/perf_cycle_codes/pmuv3.h>
//#include "/home/ubuntu/perf_cycle_codes/c_pmuv3_functions.h"

/* HEADERS FOR PMUV3 */
// SPDX-License-Identifier: GPL-2.0

extern "C"{

#include "/home/ubuntu/perf_cycle_codes/pmuv3_plugin.h"
//#include "/home/ubuntu/perf_cycle_codes/pmuv3.h"

}

/*25 MHZ timer generic */
#include <stdio.h>
#include <stdint.h>
#include <time.h> // For time measurements in C
#include <unistd.h> // For sleep functionality in C 

namespace {
std::mt19937 rgen(0);
std::string  enc_type        = "neon";
unsigned     nof_repetitions = 1000;
bool         silent          = false;
} // namespace

static void usage(const char* prog)
{
  fmt::print("Usage: {} [-R repetitions] [-s silent]\n", prog);
  fmt::print("\t-R Repetitions [Default {}]\n", nof_repetitions);
  fmt::print("\t-T Encoder type generic, avx2 or neon [Default {}]\n", enc_type);
  fmt::print("\t-s Toggle silent operation [Default {}]\n", silent);
  fmt::print("\t-h Show this message\n");
}

static void parse_args(int argc, char** argv)
{
  int opt = 0;
  while ((opt = getopt(argc, argv, "R:T:sh")) != -1) {
    switch (opt) {
      case 'R':
        nof_repetitions = std::strtol(optarg, nullptr, 10);
        break;
      case 'T':
        enc_type = std::string(optarg);
        break;
      case 's':
        silent = (!silent);
        break;
      case 'h':
      default:
        usage(argv[0]);
        exit(0);
    }
  }
}


using namespace srsran;
using namespace srsran::ldpc;
double calculatePercentile50(const std::vector<int>& data) {
    if (data.empty()) {
        // Handle the case where the data vector is empty
        return 0.0; // or some default value
    }

    // Make a copy of the data to avoid modifying the original vector
    std::vector<int> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());

    // Calculate the index corresponding to the 50th percentile
    double index = 0.5 * (sortedData.size() - 1);

    // Calculate the fractional and integer parts of the index
    double fractionalPart = std::fmod(index, 1.0);
    int integerPart = static_cast<int>(index);

    // Interpolate between the two nearest values
    int lowerValue = sortedData[integerPart];
    int upperValue = sortedData[integerPart + 1];
    double interpolatedValue = lowerValue + fractionalPart * (upperValue - lowerValue);

    return interpolatedValue;
}
double calculatePercentile(const std::vector<int>& data, double percentile) {
    if (data.empty()) {
        // Handle the case where the data vector is empty
        return 0.0; // or some default value
    }

    // Make a copy of the data to avoid modifying the original vector
    std::vector<int> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());

    // Calculate the index corresponding to the desired percentile
    double index = percentile * (sortedData.size() - 1);

    // Calculate the fractional and integer parts of the index
    double fractionalPart = std::fmod(index, 1.0);
    int integerPart = static_cast<int>(index);

    // Interpolate between the two nearest values
    int lowerValue = sortedData[integerPart];
    int upperValue = sortedData[integerPart + 1];
    double interpolatedValue = lowerValue + fractionalPart * (upperValue - lowerValue);

    return interpolatedValue;
}
//int tests_verbose;
int main(int argc, char** argv)
{
  //#ifdef __cplusplus
  //extern "C"{
  //#endif
  //__T("test evsel", !test_evsel(argc, argv));
  std::ofstream outFile("percentiles.csv");
  outFile << "Set,BG,LS,CB_LEN,50th Percentile,75th Percentile,99th Percentile,99.99th Percentile, Maximum\n";
  parse_args(argc, argv);
  
  benchmarker perf_meas_time("LDPC encoder " + enc_type + " Time", nof_repetitions); //Added by Gayathri
  benchmarker perf_meas_generic("LDPC encoder " + enc_type, nof_repetitions);
  fmt::print("\t-h Show this message TEST  Make\n");
  //uint64_t value_before = read_timer();
  //std::cout << "Value before: " << value_before << std::endl;
  // Capture the current time point
 // std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
  int iteration_counter = 0; 
  int cycle_count_prints = 0;
  std::vector<int> cd_arr;
  unsigned outst = 0, soutst = 0, inner = 0;
  std::vector<int> cd_arr1;
  std::vector<uint64_t> minValues;
  std::vector<uint64_t> maxValues;
  std::vector<double> avgValues;
  std::vector<double> percentiles50;
  std::vector<double> percentiles75;
  std::vector<double> percentiles99;
  std::vector<double> percentiles999;

  std::vector<int> bg_arr;
  std::vector<int> ls_arr;
  std::vector<int> cb_arr;
  std::vector<int> subset;
  const int subsetSize = 1000;

  for (const ldpc_base_graph_type& bg : {ldpc_base_graph_type::BG1, ldpc_base_graph_type::BG2}) {
	  outst++;
    for (const lifting_size_t& ls : all_lifting_sizes) {
      soutst++;
      std::shared_ptr<ldpc_encoder_factory> encoder_factory = create_ldpc_encoder_factory_sw(enc_type);
      TESTASSERT(encoder_factory);
      std::unique_ptr<ldpc_encoder> encoder = encoder_factory->create();
      TESTASSERT(encoder);

      // Set base-graph message and codeblock lengths.
      unsigned min_cb_length_bg = 24;
      unsigned max_cb_length_bg = 66;//66 was there in srs. i changed it to 68
  //    unsigned std_cb_length_bg = 68;
      unsigned msg_length_bg    = 22;
      if (bg == srsran::ldpc_base_graph_type::BG2) {
        min_cb_length_bg = 12;
        max_cb_length_bg = 50;
        msg_length_bg    = 10;
//	std_cb_length_bg = 52;
      }

      // Compute lifted messages and codeblock lengths.
      unsigned min_cb_length = min_cb_length_bg * ls;
      unsigned max_cb_length = max_cb_length_bg * ls;
      unsigned msg_length    = msg_length_bg * ls;

      for (unsigned cb_length : {min_cb_length, max_cb_length}) {
	inner++;
        // Generate message data.
	iteration_counter++;
        std::vector<uint8_t> data(msg_length);
        std::generate(data.begin(), data.end(), [&]() { return static_cast<uint8_t>(rgen() & 1); });

        // Generate codeblock.
        std::vector<uint8_t> codeblock(cb_length);
	//std::cout<<"cb_length="<<cb_length;
        srsran::codeblock_metadata::tb_common_metadata cfg_enc = {bg, ls};

        fmt::memory_buffer descr_buffer;
        fmt::format_to(descr_buffer, "BG={} LS={:<3} cb_len={}", bg, ls, cb_length);
	// Measure time performance
        perf_meas_time.new_measure(to_string(descr_buffer), data.size(), [&]() {
          encoder->encode(codeblock, data, cfg_enc);
          do_not_optimize(codeblock);
        }); //Gayathri
	// Measure cycle count before the encoder function
        //std::cout << "Cycle count for iteration " << iteration_counter << ": Before encoder" << std::endl;
       // uint64_t start_cycles = readCycleCount();
        //std::cout << "Cycle count for iteration " << iteration_counter << ": Before encoder" << std::endl;
	uint64_t minValue = UINT64_MAX;
	uint64_t maxValue = 0;
	double sum = 0.0;
	//uint64_t a = 0, b = 0;
	perf_meas_generic.new_measure(to_string(descr_buffer), data.size(), [&]() {
	__T("test evsel",!test_evsel(argc, argv, event_names[eventnum]));
	//__T("test evsel",!test_evsel(argc, argv));
	 uint64_t start_cycles = get_start_count(global_evsel, &global_counts); 
	 encoder->encode(codeblock, data, cfg_enc);   
	 uint64_t end_cycles = get_end_count(global_evsel, &global_counts);
	 shutdown_resources(); 
	
//	printf("Cycle count: %" PRIu64 "\n", endCount-startCount);
	//std::cout << "a " << ++a << std::endl; 
          #if 0
#endif
         // Measure cycle count after the encoder function
          uint64_t cycle_diff = end_cycles - start_cycles;

	  if (cb_length == max_cb_length){
	       //std::cout << "b " << ++b << std::endl;
	  	cd_arr.push_back(cycle_diff);
		cd_arr1.push_back(cycle_diff);
		bg_arr.push_back(static_cast<int>(bg));
		ls_arr.push_back(ls);
		cb_arr.push_back(cb_length);
		minValue = std::min(minValue, cycle_diff);
		maxValue = std::max(maxValue, cycle_diff);
		sum += cycle_diff;
	  }
	  cycle_count_prints++;
	  //std::cout << "Cycle count for iteration " << iteration_counter << ": After encoder (Print " << cycle_count_prints << ")" << std::endl;
          //std::cout << "Cycle count value: " << cycle_diff << std::endl;
	  do_not_optimize(codeblock);
	 return 0; //gayathri
        });
	double average = sum / cd_arr1.size();
	minValues.push_back(minValue);

	maxValues.push_back(maxValue);
	avgValues.push_back(average);
       // for (unsigned i = 0; i < cd_arr.size(); i++){
	//	std::cout << "Cd arr: " << cd_arr[i] << std::endl;}
	//cd_arr.clear();
        cd_arr1.clear();
	// Measure cycle count after the encoder function
       // uint64_t end_cycles = readCycleCount();
       // uint64_t cycle_diff = end_cycles - start_cycles;
        // Print cycle count and iteration information
       // cycle_count_prints++;
       // std::cout << "Cycle count for iteration " << iteration_counter << ": After encoder (Print " << cycle_count_prints << ")" << std::endl;
       // std::cout << "Cycle count value: " << cycle_diff << std::endl;
       //	std::cout << "Cycle count: " << cycle_diff << std::endl;
      }
    }
  }
  std::cout << "size of cd_arr " << cd_arr.size() << std::endl;
  for (size_t i = 0; i < cd_arr.size(); i++) {
	  subset.push_back(cd_arr[i]);
	  if (subset.size() == subsetSize) {
		//for(unsigned j =0;j<subset.size();j++){
		      //std::cout << subset[j] << " ";}
		//std::cout << std::endl;	
	  	//std::cout << "bg = " << bg_arr[i] << " ls = " << ls_arr[i] << " cb_len = " << cb_arr[i] << std::endl;
		double percentile50 = calculatePercentile(subset, 0.5);   // 50th percentile
         	double percentile75 = calculatePercentile(subset, 0.75);  // 75th percentile
         	double percentile99 = calculatePercentile(subset, 0.99);  // 99th percentile
         	double percentile9999 = calculatePercentile(subset, 0.9999); // 99.99th percentile
         	double percentile100 = calculatePercentile(subset, 1); // 100th percentile
         	std::cout << "Combination " << i / 1000 + 1 << ": 50th Percentile = " << percentile50 << std::endl;
         	std::cout << "Combination " << i / 1000 + 1 << ": 75th Percentile = " << percentile75 << std::endl;
         	std::cout << "Combination " << i / 1000 + 1 << ": 99th Percentile = " << percentile99 << std::endl;
         	std::cout << "Combination " << i / 1000 + 1 << ": 99.99th Percentile = " << percentile9999 << std::endl;
         	std::cout << "Combination " << i / 1000 + 1 << ": 100th Percentile = " << percentile100 << std::endl;
         	// Save percentiles to the CSV file
		outFile << "Set " << i << ","  << bg_arr[i] << "," << ls_arr[i] << "," << cb_arr[i] << "," << percentile50 << "," << percentile75 << "," << percentile99 << "," << percentile9999 <<  "," << percentile100 <<"\n";
		subset.clear();
}
} 

  cd_arr.clear();
  for (unsigned i = 0; i < minValues.size(); i++){
     if(i%2 !=0) 
     	std::cout << "Set " <<i<<": -> Minimum = " << minValues[i]<<", Maximum = "<<maxValues[i]<<", Average = "<<avgValues[i]<< std::endl;
  }

  perf_meas_time.print_percentiles_time("nanoseconds"); //Gayathri
  perf_meas_generic.print_percentiles_throughput("bits");
  outFile.close();

}
