#include <iostream>
#include <concepts>
#include <cstddef>
#include <optional>
#include <utility>

#include <optional>            // for std::optional
#include <xtensor/xarray.hpp>  // for ndarray

using Arr1 = xt::xarray<double, xt::layout_type::row_major>;
using Arr2 = xt::xarray<double, xt::layout_type::row_major>;

struct Options {
    size_t max_iter;
    double tol;
};

enum CutStatus {
    Success,
    NoSoln,
    NoEffect,
    SmallEnough
};

struct CInfo {
    bool feasible;
    size_t num_iters;
    CutStatus status;
};

template<typename T> using ArrayType = typename T::ArrayType;
template<typename T> using CutChoices = typename T::CutChoices;
template<typename T> using Cut = std::pair<ArrayType<T>, CutChoices<T>>;
template<typename T> using RetQ = std::tuple<Cut<T>, bool, ArrayType<T>, bool>;


template<class Oracle>
concept OracleFeas = requires(Oracle omega, const ArrayType<Oracle> &x) {
    typename Oracle::ArrayType;   // double for 1D; ndarray::Arr1 for general
    typename Oracle::CutChoices;  // double for single cut; (double, Option<double) for parallel cut
    { omega.assess_feas(x) } -> std::convertible_to<std::optional<Cut<Oracle>>>;
};

template<class Oracle>
concept OracleOptim = requires(Oracle omega, const ArrayType<Oracle> &x, double &t) {
    typename Oracle::ArrayType;   // double for 1D; ndarray::Arr1 for general
    typename Oracle::CutChoices;  // double for single cut; (double, Option<double) for parallel cut
    { omega.assess_optim(x, t) } -> std::convertible_to<std::pair<Cut<Oracle>, bool>>;
};

template<class Oracle>
concept OracleQ = requires(Oracle omega, const ArrayType<Oracle> &x, double &t, bool retry) {
    typename Oracle::ArrayType;   // double for 1D; ndarray::Arr1 for general
    typename Oracle::CutChoices;  // double for single cut; (double, Option<double) for parallel cut
    { omega.assess_q(x, t, retry) } -> std::convertible_to<RetQ<Oracle>>;
};

template<class Oracle>
concept OracleBS = requires(Oracle omega, double &t) {
    { omega.assess_bs(t) } -> std::convertible_to<bool>;
};

template<class Space, typename T>
concept SearchSpace = requires(Space ss, const std::pair<ArrayType<Space>, T> &cut) {
    typename Space::ArrayType;  // double for 1D; ndarray::Arr1 for general
    { ss.xc() } -> std::convertible_to<ArrayType<Space>>;
    { ss.update(cut) } -> std::convertible_to<std::pair<CutStatus, double>>;
};

template<typename Oracle, typename Space>
requires OracleFeas<Oracle> && SearchSpace<Space, CutChoices<Oracle>>
auto cutting_plane_feas(Oracle &omega, Space &ss, const Options &options) -> CInfo {
    for (size_t niter = 1; niter < options.max_iter; ++niter) {
        const auto cut = omega.assess_feas(ss.xc());  // query the oracle at &ss.xc()
        if (!cut) {
            // feasible sol'n obtained
            const auto[cutstatus, tsq] = ss.update(*cut);  // update ss
            if (cutstatus != CutStatus::Success) {
                return {false, niter, cutstatus};
            }
            if (tsq < options.tol) {
                return {false, niter, CutStatus::SmallEnough};
            }
        } else {
            return {true, niter, CutStatus::Success};
        }
    }
    return {false, options.max_iter, CutStatus::NoSoln};
}

template<typename T, typename Oracle>
requires OracleBS<Oracle>
auto bsearch(Oracle &omega, std::pair<T, T> &intvl, const Options &options) -> CInfo {
    auto&[lower, upper] = intvl;
    assert(lower <= upper);
    const auto u_orig = upper;

    for (size_t niter = 1; niter < options.max_iter; ++niter) {
        const auto tau = (upper - lower) / 2;  // T may be an integer
        if (tau < options.tol) {
            return {upper != u_orig, niter, CutStatus::SmallEnough};
        }
        auto t = lower;  // l may be `i32` or `Fraction`
        t += tau;
        if (omega.assess_bs(t)) {
            // feasible sol'n obtained
            upper = t;
        } else {
            lower = t;
        }
    }
    return {upper != u_orig, options.max_iter, CutStatus::NoSoln};
};

struct MyOracle {
    using ArrayType = Arr1;
    using CutChoices = double;
    using Cut = std::pair<Arr1, double>;

    /**
     * @brief
     *
     * @param[in] z
     * @return std::optional<Cut>
     */
    auto assess_feas(const Arr1 &z) -> std::optional<Cut> {
        const auto x = z[0];
        const auto y = z[1];

        // constraint 1: x + y <= 3
        const auto fj = x + y - 3.0;
        if (fj > 0.0) {
            return {{Arr1{1.0, 1.0}, fj}};
        }
        // constraint 2: x - y >= 1
        const auto fj2 = -x + y + 1.0;
        if (fj2 > 0.0) {
            return {{Arr1{-1.0, 1.0}, fj2}};
        }
        return {};
    }
};

struct MyOracleBS {
    double f;
    auto assess_bs(double &z) -> bool {
        auto res = z<f;
        f = z;
        return res;
    }
};

class EllCalc {
private:
    double n_float;
    double n_plus_1;
    double half_n;
    double n_sq;
    double c1;
    double c2;
    double c3;

public:
    double rho;
    double sigma;
    double delta;
    double tsq;
    bool use_parallel_cut;

    /**
     * @brief Construct a new Ell Calc object
     *
     * @param n_float
     */
    EllCalc(double n_float)
            : n_plus_1{n_float + 1.0},
              half_n{n_float / 2.0},
              n_sq{n_float * n_float},
              c1{n_sq / (n_sq - 1.0)},
              c2{2.0 / n_plus_1},
              c3{n_float / n_plus_1},
              rho{0.0},
              sigma{0.0},
              delta{0.0},
              tsq{0.0},
              use_parallel_cut{true} {}

    /**
     * @brief
     *
     * @param[in] b0
     * @param[in] b1
     * @return CutStatus
     */
    auto calc_ll_core(double b0, double b1) -> CutStatus { return CutStatus::Success; }

    /**
     * @brief
     *
     * @param[in] b1
     * @param[in] b1sq
     */
    auto calc_ll_cc(double b1, double b1sqn) -> void {};

    /**
     * @brief Deep Cut
     *
     * @param[in] beta
     * @return CutStatus
     */
    auto calc_dc(double beta) -> CutStatus { return CutStatus::Success; }

    /**
     * @brief Central Cut
     *
     * @param[in] tau
     */
    auto calc_cc(double tau) -> void {};

    /**
     * @brief Get the results object
     *
     * @return std::array<double, 4>
     */
    auto get_results() const -> std::array<double, 4> {
        return {this->rho, this->sigma, this->delta, this->tsq};
    }
};

class Ell {
    using Self = Ell;
    using Parallel = std::pair<double, std::optional<double>>;

    size_t n;
    double kappa;
    Arr2 mq;
    Arr1 xc_;
    EllCalc helper;

public:
    using ArrayType = Arr1;

    bool no_defer_trick;

    /**
     * @brief Construct a new Ell Stable object
     *
     * @param[in] kappa
     * @param[in] mq
     * @param[in] xc
     */
    Ell(double kappa, Arr2 mq, Arr1 xc)
            : n{xc.size()},
              kappa{kappa},
              mq{std::move(mq)},
              xc_{std::move(xc)},
              helper(double(n)),
              no_defer_trick{false} {}

    /**
     * @brief Construct a new Ell Stable object
     *
     * @param[in] val
     * @param[in] xc
     */
    Ell(Arr1 val, Arr1 xc) : Ell{1.0, xt::diag(val), std::move(xc)} {}

    /**
     * @brief Construct a new Ell Stable object
     *
     * @param[in] alpha
     * @param[in] xc
     */
    Ell(double alpha, Arr1 xc) : Ell{alpha, xt::eye(xc.size()), std::move(xc)} {}

    /**
     * @brief
     *
     * @param[in] grad
     * @param[in] beta
     * @return std::pair<CutStatus, double>
     */
    auto update_single(const Arr1 &grad, const double &beta) -> std::pair<CutStatus, double> {
        // const auto [grad, beta] = cut;
        auto mq_g = Arr1{xt::zeros<double>({this->n})};  // initial x0
        auto omega = 0.0;
        for (auto i = 0; i < this->n; i++) {
            for (auto j = 0; j < this->n; j++) {
                mq_g[i] += this->mq[{i, j}] * grad[j];
            }
            omega += mq_g[i] * grad[i];
        }

        this->helper.tsq = this->kappa * omega;
        const auto status = this->helper.calc_dc(beta);
        if (status != CutStatus::Success) {
            return {status, this->helper.tsq};
        }

        this->xc_ -= (this->helper.rho / omega) * mq_g;  // n

        const auto r = this->helper.sigma / omega;
        for (auto i = 0; i < this->n; i++) {
            const auto r_mq_g = r * mq_g[i];
            for (auto j = 0; j < i; j++) {
                this->mq[{i, j}] -= r_mq_g * mq_g[j];
                this->mq[{j, i}] = this->mq[{i, j}];
            }
            this->mq[{i, i}] -= r_mq_g * mq_g[i];
        }

        this->kappa *= this->helper.delta;

        if (this->no_defer_trick) {
            this->mq *= this->kappa;
            this->kappa = 1.0;
        }
        return {status, this->helper.tsq};
    }

/**
 * @brief Update ellipsoid core function using the cut
 *
 *  $grad^T * (x - xc) + beta <= 0$
 *
 * @tparam T
 * @param[in] cut
 * @return (i32, double)
 */
    auto update_parallel(const Arr1 &grad, const std::pair<double, std::optional<double>> &beta)
    -> std::pair<CutStatus, double> {
        // const auto [grad, beta] = cut;
        auto mq_g = Arr1{xt::zeros<double>({this->n})};  // initial x0
        auto omega = 0.0;
        for (auto i = 0; i < this->n; i++) {
            for (auto j = 0; j < this->n; j++) {
                mq_g[i] += this->mq[{i, j}] * grad[j];
            }
            omega += mq_g[i] * grad[i];
        }

        this->helper.tsq = this->kappa * omega;
        const auto[b0, b1_opt] = beta;
        const auto status = b1_opt ? this->helper.calc_ll_core(b0, *b1_opt) : this->helper.calc_dc(b0);
        if (status != CutStatus::Success) {
            return {status, this->helper.tsq};
        }

        this->xc_ -= (this->helper.rho / omega) * mq_g;  // n

        // n*(n+1)/2 + n
        // this->mq -= (this->sigma / omega) * xt::linalg::outer(mq_g, mq_g);

        const auto r = this->helper.sigma / omega;
        for (auto i = 0; i < this->n; i++) {
            const auto r_mq_g = r * mq_g[i];
            for (auto j = 0; j < this->n; j++) {
                this->mq[{i, j}] -= r_mq_g * mq_g[j];
                this->mq[{j, i}] = this->mq[{i, j}];
            }
            this->mq[{i, i}] -= r_mq_g * mq_g[i];
        }

        this->kappa *= this->helper.delta;

        if (this->no_defer_trick) {
            this->mq *= this->kappa;
            this->kappa = 1.0;
        }
        return {status, this->helper.tsq};
    }

    /**
     * @brief copy the whole array anyway
     *
     * @return Arr1
     */
    auto xc() const -> Self::ArrayType { return this->xc_; }

    /**
     * @brief
     *
     * @tparam T
     * @param[in] cut
     * @return std::pair<CutStatus, double>
     */
    template<typename T>
    auto update(const std::pair<Self::ArrayType, T> &cut)
    -> std::pair<CutStatus, double> {
        const auto[grad, beta] = cut;
        if constexpr (std::is_same_v<T, double>) {
            return this->update_single(grad, beta);
        } else if constexpr (std::is_same_v<T, Parallel>) {
            return this->update_parallel(grad, beta);
        } else {
            // static_assert(false, "Not supported type");
            return {CutStatus::NoSoln, 0.0};
        }
    }
};

int main() {
    auto ell = Ell(Arr1{10.0, 10.0}, Arr1{0.0, 0.0});
    auto oracle = MyOracle{};
    const auto options = Options{2000, 1e-12};
    const auto[feasible, _niter, _status] = cutting_plane_feas(oracle, ell, options);
    std::cout << feasible << _niter << _status << std::endl;
    auto pair = std::make_pair(1e1,2e1);
    auto oracle_ = MyOracleBS{};
    const auto info = bsearch(oracle_,pair , options);
    std::cout << feasible << std::endl;
    return 0;
}
