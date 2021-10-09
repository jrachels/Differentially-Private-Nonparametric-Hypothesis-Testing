// @author Jeremy Rachels
// Contact: jeremyrachels@gmail.com

#ifndef DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEWILCOXON_H
#define DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEWILCOXON_H

// TO DO: make work with any input that is totally ordered and has a minus operation

#include <vector>
#include <utility>
#include <cassert>
#include <map>
#include <random>
#include <concepts>
#include <cmath>

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/laplacedistribution/laplacedistribution.h"

namespace dpnht {

	template<std::uniform_random_bit_generator Engine = std::mt19937>
	[[nodiscard]] std::pair<const double, const double> PrivateWilcoxon(const std::vector<std::pair<double, double>>& dataset, const double epsilon, Engine& Eng, const int p_value_trials);


	[[nodiscard]] std::pair<const double, const double> PrivateWilcoxon(const std::vector<std::pair<double, double>>& dataset, const double epsilon, const int p_value_trials);

	namespace internal_wilcoxon {

		// Pratt Variant
		double PublicWilcoxon(const std::vector<std::pair<double, double>>& dataset);

		template<std::uniform_random_bit_generator Engine>
		double PrivateWilcoxonStatistic(const std::vector<std::pair<double, double>>& dataset, double epsilon, Engine& Eng);


		template<std::uniform_random_bit_generator Engine>
		double ComputePValue(const int dataset_size, const double private_w_statistic, const double epsilon, Engine& Eng, const int p_value_trials);

		template<std::uniform_random_bit_generator Engine = std::mt19937>
		std::pair<const double, const double> PrivateWilcoxonHelper(const std::vector<std::pair<double, double>>& dataset, const double epsilon, Engine& Eng, const int p_value_trials);

	}
}

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/privatewilcoxon/privatewilcoxon.tpp"

#endif