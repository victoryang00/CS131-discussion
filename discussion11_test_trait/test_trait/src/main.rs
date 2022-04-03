use ndarray::prelude::*;
type CInfo = (bool, usize, CutStatus);

pub struct Options {
    pub max_iter: usize,
    pub tol: f64,
}


#[derive(Debug, PartialEq, Eq)]
pub enum CutStatus {
    Success,
    NoSoln,
    NoEffect,
    SmallEnough,
}

pub trait UpdateByCutChoices<SS> {
    type ArrayType; // f64 for 1D; ndarray::Array1<f64> for general

    fn update_by(&self, ss: &mut SS, grad: &Self::ArrayType) -> (CutStatus, f64);
}

/// Oracle for feasibility problems
pub trait OracleFeas {
    type ArrayType;
    // f64 for 1D; ndarray::Array1<f64> for general
    type CutChoices;
    // f64 for single cut; (f64, Option<f64) for parallel cut
    fn assess_feas(&mut self, x: &Self::ArrayType) -> Option<(Self::ArrayType, Self::CutChoices)>;
}

/// Oracle for optimization problems
pub trait OracleOptim {
    type ArrayType;
    // f64 for 1D; ndarray::Array1<f64> for general
    type CutChoices;
    // f64 for single cut; (f64, Option<f64) for parallel cut
    fn assess_optim(
        &mut self,
        x: &Self::ArrayType,
        t: &mut f64,
    ) -> ((Self::ArrayType, Self::CutChoices), bool);
}

/// Oracle for quantized optimization problems
pub trait OracleQ {
    type ArrayType;
    // f64 for 1D; ndarray::Array1<f64> for general
    type CutChoices;
    // f64 for single cut; (f64, Option<f64) for parallel cut
    fn assess_q(
        &mut self,
        x: &Self::ArrayType,
        t: &mut f64,
        retry: bool,
    ) -> (
        (Self::ArrayType, Self::CutChoices),
        bool,
        Self::ArrayType,
        bool,
    );
}

pub trait OracleBS {
    fn assess_bs(&mut self, t: f64) -> bool;
}


pub trait SearchSpace {
    type ArrayType;
    // f64 for 1D; ndarray::Array1<f64> for general
    fn xc(&self) -> Self::ArrayType;
    fn update<T>(&mut self, cut: &(Self::ArrayType, T)) -> (CutStatus, f64)
        where
            T: UpdateByCutChoices<Self, ArrayType=Self::ArrayType>,
            Self: Sized;
}

#[allow(dead_code)]
pub fn cutting_plane_feas<T, Oracle, Space>(
    omega: &mut Oracle,
    ss: &mut Space,
    options: &Options,
) -> CInfo
    where
        T: UpdateByCutChoices<Space, ArrayType=Oracle::ArrayType>,
        Oracle: OracleFeas<CutChoices=T>,
        Space: SearchSpace<ArrayType=Oracle::ArrayType>,
{
    for niter in 1..options.max_iter {
        let cut_option = omega.assess_feas(&ss.xc()); // query the oracle at &ss.xc()
        if let Some(cut) = cut_option {
            // feasible sol'n obtained
            let (cutstatus, tsq) = ss.update::<T>(&cut); // update ss
            if cutstatus != CutStatus::Success {
                return (false, niter, cutstatus);
            }
            if tsq < options.tol {
                return (false, niter, CutStatus::SmallEnough);
            }
        } else {
            return (true, niter, CutStatus::Success);
        }
    }
    (false, options.max_iter, CutStatus::NoSoln)
}

#[allow(dead_code)]
pub fn cutting_plane_optim<T, Oracle, Space>(
    omega: &mut Oracle,
    ss: &mut Space,
    t: &mut f64,
    options: &Options,
) -> (Option<Oracle::ArrayType>, usize, CutStatus)
    where
        T: UpdateByCutChoices<Space, ArrayType=Oracle::ArrayType>,
        Oracle: OracleOptim<CutChoices=T>,
        Space: SearchSpace<ArrayType=Oracle::ArrayType>,
{
    let mut x_best: Option<Oracle::ArrayType> = None;
    let mut status = CutStatus::NoSoln;

    for niter in 1..options.max_iter {
        let (cut, shrunk) = omega.assess_optim(&ss.xc(), t); // query the oracle at &ss.xc()
        if shrunk {
            // best t obtained
            x_best = Some(ss.xc());
            status = CutStatus::Success;
        }
        let (cutstatus, tsq) = ss.update::<T>(&cut); // update ss
        if cutstatus != CutStatus::Success {
            return (x_best, niter, cutstatus);
        }
        if tsq < options.tol {
            return (x_best, niter, CutStatus::SmallEnough);
        }
    }
    (x_best, options.max_iter, status)
}

#[allow(dead_code)]
pub fn cutting_plane_q<T, Oracle, Space>(
    omega: &mut Oracle,
    ss: &mut Space,
    t: &mut f64,
    options: &Options,
) -> (Option<Oracle::ArrayType>, usize, CutStatus)
    where
        T: UpdateByCutChoices<Space, ArrayType=Oracle::ArrayType>,
        Oracle: OracleQ<CutChoices=T>,
        Space: SearchSpace<ArrayType=Oracle::ArrayType>,
{
    let mut x_best: Option<Oracle::ArrayType> = None;
    let mut status = CutStatus::NoSoln; // note!!!
    let mut retry = false;

    for niter in 1..options.max_iter {
        let (cut, shrunk, x0, more_alt) = omega.assess_q(&ss.xc(), t, retry); // query the oracle at &ss.xc()
        if shrunk {
            // best t obtained
            x_best = Some(x0); // x0
        }
        let (cutstatus, tsq) = ss.update::<T>(&cut); // update ss
        match &cutstatus {
            CutStatus::NoEffect => {
                if !more_alt {
                    // more alt?
                    return (x_best, niter, status);
                }
                status = cutstatus;
                retry = true;
            }
            CutStatus::NoSoln => {
                return (x_best, niter, CutStatus::NoSoln);
            }
            _ => {}
        }
        if tsq < options.tol {
            return (x_best, niter, CutStatus::SmallEnough);
        }
    }
    (x_best, options.max_iter, status)
}

#[allow(dead_code)]
pub fn bsearch<Oracle>(omega: &mut Oracle, intvl: &mut (f64, f64), options: &Options) -> CInfo
    where
        Oracle: OracleBS,
{
    // assume monotone
    // auto& [lower, upper] = I;
    let &mut (mut lower, mut upper) = intvl;
    assert!(lower <= upper);
    let u_orig = upper;
    let mut status = CutStatus::NoSoln;

    let mut niter = 0;
    while niter < options.max_iter {
        niter += 1;
        let tau = (upper - lower) / 2.0;
        if tau < options.tol {
            status = CutStatus::SmallEnough;
            break;
        }
        let mut t = lower; // l may be `i32` or `Fraction`
        t += tau;
        if omega.assess_bs(t) {
            // feasible sol'n obtained
            upper = t;
        } else {
            lower = t;
        }
    }
    (upper != u_orig, niter, status)
}

struct MyOracleBS {
    d: f64,
}

impl OracleBS for MyOracleBS {
    fn assess_bs(&mut self, t: f64) -> bool {
        let res = t < self.d;
        self.d = t;
        res
    }
}

struct MyOracle {}

impl OracleFeas for MyOracle {
    type ArrayType =  Array1<f64>;
    type CutChoices = f64;
    /**
     * @brief
     *
     * @param[in] z
     * @return std::optional<Cut>
     */
      fn assess_feas(&mut self, z: &Array1<f64>) -> Option<(Array1<f64>, f64)> {
        let x = z[0];
        let y = z[1];

        // constraint 1: x + y <= 3
        let fj = x + y - 3.0;
        if fj > 0.0 {
            return Some((array![1.0, 1.0], fj));
        }
        // constraint 2: x - y >= 1
        let fj = -x + y + 1.0;
        if fj > 0.0 {
            return Some((array![-1.0, 1.0], fj));
        }
        return None;
    }
}

/**
 * @brief Ellipsoid Search Space
 *
 *  Ell = {x | (x - xc)^T mq^-1 (x - xc) \le \kappa}
 *
 * Keep $mq$ symmetric but no promise of positive definite
 */
#[derive(Debug, Clone)]
pub struct Ell {
    pub no_defer_trick: bool,

    mq: Array2<f64>,
    xc: Array1<f64>,
    kappa: f64,
    n: usize,
    helper: EllCalc,
}

impl Ell {
    /**
     * @brief Construct a new Ell object
     *
     * @tparam V
     * @tparam U
     * @param kappa
     * @param mq
     * @param x
     */
    pub fn new_with_matrix(kappa: f64, mq: Array2<f64>, xc: Array1<f64>) -> Ell {
        let n = xc.len();
        let helper = EllCalc::new(n as f64);

        Ell {
            kappa,
            mq,
            xc,
            n,
            helper,
            no_defer_trick: false,
        }
    }

    /**
     * @brief Construct a new Ell object
     *
     * @param[in] val
     * @param[in] x
     */
    pub fn new(val: Array1<f64>, xc: Array1<f64>) -> Ell {
        Ell::new_with_matrix(1.0, Array2::from_diag(&val), xc)
    }

    // /**
    //  * @brief Set the xc object
    //  *
    //  * @param[in] xc
    //  */
    // pub fn set_xc(&mut self, xc: Arr) { self.xc = xc; }

    /**
     * @brief Update ellipsoid core function using the cut
     *
     *  $grad^T * (x - xc) + beta <= 0$
     *
     * @tparam T
     * @param[in] cut
     * @return (i32, f64)
     */
    fn update_single(&mut self, grad: &Array1<f64>, beta: &f64) -> (CutStatus, f64) {
        // let (grad, beta) = cut;
        let mut mq_g = Array1::zeros(self.n); // initial x0
        let mut omega = 0.0;
        for i in 0..self.n {
            for j in 0..self.n {
                mq_g[i] += self.mq[[i, j]] * grad[j];
            }
            omega += mq_g[i] * grad[i];
        }

        self.helper.tsq = self.kappa * omega;
        let status = self.helper.calc_dc(*beta);
        if status != CutStatus::Success {
            return (status, self.helper.tsq);
        }

        self.xc -= &((self.helper.rho / omega) * &mq_g); // n

        // n*(n+1)/2 + n
        // self.mq -= (self.sigma / omega) * xt::linalg::outer(mq_g, mq_g);

        let r = self.helper.sigma / omega;
        for i in 0..self.n {
            let r_mq_g = r * mq_g[i];
            for j in 0..i {
                self.mq[[i, j]] -= r_mq_g * mq_g[j];
                self.mq[[j, i]] = self.mq[[i, j]];
            }
            self.mq[[i, i]] -= r_mq_g * mq_g[i];
        }

        self.kappa *= self.helper.delta;

        if self.no_defer_trick {
            self.mq *= self.kappa;
            self.kappa = 1.0;
        }
        (status, self.helper.tsq)
    }

    /**
     * @brief Update ellipsoid core function using the cut
     *
     *  $grad^T * (x - xc) + beta <= 0$
     *
     * @tparam T
     * @param[in] cut
     * @return (i32, f64)
     */
    fn update_parallel(
        &mut self,
        grad: &Array1<f64>,
        beta: &(f64, Option<f64>),
    ) -> (CutStatus, f64) {
        // let (grad, beta) = cut;
        let mut mq_g = Array1::zeros(self.n); // initial x0
        let mut omega = 0.0;
        for i in 0..self.n {
            for j in 0..self.n {
                mq_g[i] += self.mq[[i, j]] * grad[j];
            }
            omega += mq_g[i] * grad[i];
        }

        self.helper.tsq = self.kappa * omega;

        let (b0, b1_opt) = *beta;
        let status = if let Some(b1) = b1_opt {
            self.helper.calc_ll_core(b0, b1)
        } else {
            self.helper.calc_dc(b0)
        };
        if status != CutStatus::Success {
            return (status, self.helper.tsq);
        }

        self.xc -= &((self.helper.rho / omega) * &mq_g); // n

        // n*(n+1)/2 + n
        // self.mq -= (self.sigma / omega) * xt::linalg::outer(mq_g, mq_g);

        let r = self.helper.sigma / omega;
        for i in 0..self.n {
            let r_mq_g = r * mq_g[i];
            for j in 0..i {
                self.mq[[i, j]] -= r_mq_g * mq_g[j];
                self.mq[[j, i]] = self.mq[[i, j]];
            }
            self.mq[[i, i]] -= r_mq_g * mq_g[i];
        }

        self.kappa *= self.helper.delta;

        if self.no_defer_trick {
            self.mq *= self.kappa;
            self.kappa = 1.0;
        }
        (status, self.helper.tsq)
    }
}

impl SearchSpace for Ell {
    type ArrayType = Array1<f64>;

    /**
     * @brief copy the whole array anyway
     *
     * @return Arr
     */
    fn xc(&self) -> Self::ArrayType {
        self.xc.clone()
    }

    fn update<T>(&mut self, cut: &(Self::ArrayType, T)) -> (CutStatus, f64)
        where
            T: UpdateByCutChoices<Self, ArrayType = Self::ArrayType>,
    {
        let (grad, beta) = cut;
        beta.update_by(self, grad)
    }
}

impl UpdateByCutChoices<Ell> for f64 {
    type ArrayType = Array1<f64>;

    fn update_by(&self, ell: &mut Ell, grad: &Self::ArrayType) -> (CutStatus, f64) {
        let beta = self;
        ell.update_single(grad, beta)
    }
}

impl UpdateByCutChoices<Ell> for (f64, Option<f64>) {
    type ArrayType = Array1<f64>;

    fn update_by(&self, ell: &mut Ell, grad: &Self::ArrayType) -> (CutStatus, f64) {
        let beta = self;
        ell.update_parallel(grad, beta)
    }
}

/**
 * @brief Ellipsoid Search Space
 *
 *  EllCalc = {x | (x - xc)^T mq^-1 (x - xc) \le \kappa}
 *
 * Keep $mq$ symmetric but no promise of positive definite
 */
#[derive(Debug, Clone)]
pub struct EllCalc {
    pub use_parallel_cut: bool,

    pub rho: f64,
    pub sigma: f64,
    pub delta: f64,
    pub tsq: f64,

    n_float: f64,
    n_plus_1: f64,
    half_n: f64,
    c1: f64,
    c2: f64,
    c3: f64,
}

impl EllCalc {
    /**
     * @brief Construct a new EllCalc object
     *
     * @tparam V
     * @tparam U
     * @param kappa
     * @param mq
     * @param x
     */
    pub fn new(n_float: f64) -> EllCalc {
        let n_plus_1 = n_float + 1.0;
        let half_n = n_float / 2.0;
        let n_sq = n_float * n_float;
        let c1 = n_sq / (n_sq - 1.0);
        let c2 = 2.0 / n_plus_1;
        let c3 = n_float / n_plus_1;

        EllCalc {
            n_float,
            n_plus_1,
            half_n,
            c1,
            c2,
            c3,
            rho: 0.0,
            sigma: 0.0,
            delta: 0.0,
            tsq: 0.0,
            use_parallel_cut: true,
        }
    }

    // pub fn update_cut(&mut self, beta: f64) -> CutStatus { self.calc_dc(beta) }

    /**
     * @brief
     *
     * @param[in] b0
     * @param[in] b1
     * @return i32
     */
    pub fn calc_ll_core(&mut self, b0: f64, b1: f64) -> CutStatus {
        // let b1sq = b1 * b1;
        let b1sqn = b1 * (b1 / self.tsq);
        let t1n = 1.0 - b1sqn;
        if t1n < 0.0 || !self.use_parallel_cut {
            return self.calc_dc(b0);
        }

        let bdiff = b1 - b0;
        if bdiff < 0.0 {
            return CutStatus::NoSoln; // no sol'n
        }

        if b0 == 0.0 {
            // central cut
            self.calc_ll_cc(b1, b1sqn);
            return CutStatus::Success;
        }

        let b0b1n = b0 * (b1 / self.tsq);
        if self.n_float * b0b1n < -1.0 {
            return CutStatus::NoEffect; // no effect
        }

        // let t0 = self.tsq - b0 * b0;
        let t0n = 1.0 - b0 * (b0 / self.tsq);
        // let t1 = self.tsq - b1sq;
        let bsum = b0 + b1;
        let bsumn = bsum / self.tsq;
        let bav = bsum / 2.0;
        let tempn = self.half_n * bsumn * bdiff;
        let xi = (t0n * t1n + tempn * tempn).sqrt();
        self.sigma = self.c3 + (1.0 - b0b1n - xi) / (bsumn * bav) / self.n_plus_1;
        self.rho = self.sigma * bav;
        self.delta = self.c1 * ((t0n + t1n) / 2.0 + xi / self.n_float);
        CutStatus::Success
    }

    /**
     * @brief
     *
     * @param[in] b1
     * @param[in] b1sq
     * @return void
     */
    pub fn calc_ll_cc(&mut self, b1: f64, b1sqn: f64) {
        let temp = self.half_n * b1sqn;
        let xi = (1.0 - b1sqn + temp * temp).sqrt();
        self.sigma = self.c3 + self.c2 * (1.0 - xi) / b1sqn;
        self.rho = self.sigma * b1 / 2.0;
        self.delta = self.c1 * (1.0 - b1sqn / 2.0 + xi / self.n_float);
    }

    /**
     * @brief Deep Cut
     *
     * @param[in] beta
     * @return i32
     */
    pub fn calc_dc(&mut self, beta: f64) -> CutStatus {
        let tau = (self.tsq).sqrt();

        let bdiff = tau - beta;
        if bdiff < 0.0 {
            return CutStatus::NoSoln; // no sol'n
        }

        if beta == 0.0 {
            self.calc_cc(tau);
            return CutStatus::Success;
        }

        let gamma = tau + self.n_float * beta;
        if gamma < 0.0 {
            return CutStatus::NoEffect; // no effect
        }

        // self.mu = (bdiff / gamma) * self.half_n_minus_1;
        self.rho = gamma / self.n_plus_1;
        self.sigma = 2.0 * self.rho / (tau + beta);
        self.delta = self.c1 * (1.0 - beta * (beta / self.tsq));
        CutStatus::Success
    }

    /**
     * @brief Central Cut
     *
     * @param[in] tau
     * @return i32
     */
    pub fn calc_cc(&mut self, tau: f64) {
        // self.mu = self.half_n_minus_1;
        self.sigma = self.c2;
        self.rho = tau / self.n_plus_1;
        self.delta = self.c1;
    }

    pub fn get_results(&self) -> [f64; 4] {
        [self.rho, self.sigma, self.delta, self.tsq]
    }
}

fn main() {
    let mut options = Options{
        max_iter: 100,
        tol: 1e-6,
    };

    let mut ell = Ell::new(array![10.0, 10.0], array![0.0, 0.0]);
    let mut omega = MyOracle {};
    let (x, niter, status) = cutting_plane_feas(&mut omega, &mut ell,  &options);
    println!("{:?}", x);
    println!("{:?}", niter);
    println!("{:?}", status);
    let mut omega_ = MyOracleBS {d:1.0};
    let (b, niter, status) = bsearch(&mut omega_, &mut (1.0, 2.0), &options);
    println!("{:?}", b);
}
