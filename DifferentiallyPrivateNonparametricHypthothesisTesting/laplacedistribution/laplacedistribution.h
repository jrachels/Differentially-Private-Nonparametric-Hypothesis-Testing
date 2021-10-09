// @author Jeremy Rachels
// Contact: jeremyrachels@gmail.com

#ifndef DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_LAPLACE_DISTRIBUTION_H
#define DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_LAPLACE_DISTRIBUTION_H

#include <concepts>
#include <random>
#include <limits>
#include <cmath>

namespace dpnht {


    // CLASS TEMPLATE laplace_distribution
    class laplace_distribution { // laplace distribution
    public:
        struct param_type { // parameter package
            param_type() {
                _Init(0.0, 1.0);
            }

            explicit param_type(double _Mean0, double _Scale0 = 1.0) {
                _Init(_Mean0, _Scale0);
            }

            _NODISCARD bool operator==(const param_type& _Right) const {
                return _Mean == _Right._Mean && _Scale == _Right._Scale;
            }

            _NODISCARD bool operator!=(const param_type& _Right) const {
                return !(*this == _Right);
            }

            _NODISCARD double mean() const {
                return _Mean;
            }

            _NODISCARD double scale() const {
                return _Scale;
            }

            void _Init(double _Mean0, double _Scale0) { // set internal state
                _STL_ASSERT(0.0 <= _Scale0, "invalid scale argument for laplace_distribution");
                _Mean = _Mean0;
                _Scale = _Scale0;
            }

            double _Mean;
            double _Scale;
        };

        laplace_distribution() : _Par(0.0, 1.0) {}

        explicit laplace_distribution(double _Mean0, double _Scale0 = 1.0) : _Par(_Mean0, _Scale0) {}

        explicit laplace_distribution(const param_type& _Par0) : _Par(_Par0) {}

        _NODISCARD double mean() const {
            return _Par.mean();
        }

        _NODISCARD double sigma() const {
            return _Par.scale();
        }

        _NODISCARD param_type param() const {
            return _Par;
        }

        void param(const param_type& _Par0) { // set parameter package
            _Par = _Par0;
        }

        _NODISCARD double(min)() const { // get smallest possible result
            return std::numeric_limits<double>::denorm_min();
        }

        _NODISCARD double(max)() const { // get largest possible result
            return std::numeric_limits<double>::max();
        }

        template <std::uniform_random_bit_generator _Engine>
        _NODISCARD double operator()(_Engine& _Eng) {
            return _Eval(_Eng, _Par);
        }

        template <std::uniform_random_bit_generator _Engine>
        _NODISCARD double operator()(_Engine& _Eng, const param_type& _Par0) {
            return _Eval(_Eng, _Par0);
        }

    private:
        template <std::uniform_random_bit_generator _Engine>
        double _Eval(_Engine& _Eng, const param_type& _Par0) { // compute next value
            double _Res;

            double random_value = std::generate_canonical<double, static_cast<size_t>(std::numeric_limits<double>::digits)>(_Eng);
            if (random_value < .5) {
                _Res = std::log(2.0 * random_value);
            }
            else {
                _Res = -std::log(2.0 - 2.0 * random_value);
            }

            return _Res * _Par0._Scale + _Par0._Mean;
        }

        param_type _Par;
    };


}

#endif