// @author Jeremy Rachels
// Contact: jeremyrachels@gmail.com

#ifndef DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEKRUSKALWALLIS_H
#define DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEKRUSKALWALLIS_H

// TO DO: can eliminate create of rankings_data object and instead work directly with the std::map from. needs to be tested. compromises clarity of code
// TO DO: instead of allowing ties, randomly break ties. will require figuring out theory behind how this decision affects p-value
// TO DO: in rankings_data in PublicKruskalWallis, for each group, find the point where the ranking exceeds the mean rank for the whole data set.
// this allows you to not need to use absolute value. 
// TO DO: move all of the static functions to cpp file
// TO DO: don't recompute values like group sizes. pass those values through function arguments.
// TO DO: replace map in computeRankings with a multi-map
// TO DO: Change warning to issue a compiler warning instead of using std::cout
// TO DO: The two PrivateKruskalWallis functions are idential except that one creates an Eng if it wasn't passed as an argument.
// create a PrivateKruskalWallisHelper function for this duplicate code.

#include <vector>
#include <cmath> 
#include <random>
#include <utility>

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/laplacedistribution/laplacedistribution.h"

namespace dpnht {

	template<std::totally_ordered datatype = double, std::uniform_random_bit_generator Engine = std::mt19937>
	[[nodiscard]] std::pair<const double, const double> PrivateKruskalWallis(const std::vector<std::vector<datatype>>& raw_data, const double epsilon, Engine& Eng, const int p_value_trials, bool group_sizes_known = false);


	template<std::totally_ordered datatype = double>
	[[nodiscard]] std::pair<const double, const double> PrivateKruskalWallis(const std::vector<std::vector<datatype>>& raw_data, const double epsilon, const int p_value_trials, bool group_sizes_known = false);


	// helper functions below
	namespace internal_kruskal_wallis {

		double KruskalWallisAbs(const std::vector<std::vector<double>>& rankings_data);


		template<std::uniform_random_bit_generator Engine>
		double PrivateKruskalWallisStatistic(const std::vector<std::vector<double>>& rankings_data, double epsilon, Engine& Eng);


		template<std::uniform_random_bit_generator Engine>
		std::vector<std::vector<double>> BuildRandomDatabase(const std::size_t num_groups, const std::size_t num_data_points, Engine& Eng);



		template<std::uniform_random_bit_generator Engine>
		std::vector<std::vector<double>> BuildRandomDatabase(const std::vector<int>& group_sizes, Engine& Eng);


		template<std::totally_ordered datatype = double>
		std::vector<std::vector<double>> ComputeRankings(const std::vector<std::vector<datatype>>& raw_data);



		// do a different p-value statistic if group sizes are known
		template<std::uniform_random_bit_generator Engine>
		double ComputePValue(const std::vector<int>& group_sizes, const double private_h_statistic, const double epsilon, Engine& Eng, const int p_value_trials);

		// Note: this provides a conservative estimate of the p-value without revealing group sizes
		template<std::uniform_random_bit_generator Engine>
		double ComputePValue(const int num_groups, const int num_data_points, const double private_h_statistic, const double epsilon, Engine& Eng, const int p_value_trials);

	}

}

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/privatekruskalwallis/privatekruskalwallis.tpp"

#endif