/*
 * problem.cpp
 *
 *  Created on: May 29, 2015
 *      Author: dmarce1
 */

#include "grhd.hpp"

#include "defs.hpp"
#include "options.hpp"
#include "grid.hpp"
#include <cmath>
#include <utility>
#include <array>

static const auto eta = simd_vector(0.5);
static const auto kap1 = simd_vector(0.5);
static const auto kap2 = simd_vector(-0.5);

constexpr integer grhd::NF;
static constexpr integer D_i = 0;
static constexpr integer tau_i = 1;
static constexpr integer Sx_i = 2;
static constexpr integer Sy_i = 3;
static constexpr integer Sz_i = 4;

static constexpr integer rho_i = 0;
static constexpr integer eps_i = 1;
static constexpr integer ux_i = 2;
static constexpr integer uy_i = 3;
static constexpr integer uz_i = 4;

static const simd_vector zero(0.0);
static const simd_vector one(1.0);
static const simd_vector two(2.0);
static const simd_vector half(0.5);
static const simd_vector fgamma(5.0 / 3.0);
static const simd_vector mu = two;

namespace gr {

constexpr real delta(integer i, integer j) {
	return i == j ? 1.0 : 0.0;
}

using scalar = simd_vector;

class vector {
private:
	static constexpr integer NDIM = 3;
	std::array<simd_vector, NDIM> data;
public:
	vector() = default;
	vector(const vector&) = default;
	vector(vector&&) = default;
	~vector() = default;
	vector& operator=(const vector&) = default;
	vector& operator=(vector &&) = default;
	const simd_vector& operator()(integer i) const {
		return data[i];
	}
	simd_vector& operator()(integer i) {
		return data[i];
	}
	vector& operator*=(const scalar& other) {
		for (integer i = 0; i != NDIM; ++i) {
			data[i] *= other;
		}
		return *this;
	}
	vector& operator/=(const scalar& other) {
		const auto other_inv = simd_vector(1) / other;
		return operator*=(other_inv);
	}
	vector operator*(const scalar& other) const {
		vector c(*this);
		for (integer i = 0; i != NDIM; ++i) {
			c.data[i] *= other;
		}
		return c;
	}
	vector operator/(const scalar& other) const {
		const auto other_inv = simd_vector(1) / other;
		return operator*(other_inv);
	}
	scalar operator*(const vector& other) const {
		scalar sum(data[0] * other(0));
		sum += data[1] * other(1);
		sum += data[2] * other(2);
		return sum;
	}
};

scalar abs(const vector& A) {
	return sqrt(A * A);
}

class tensor {
private:
	constexpr static integer SIZE = 6;
	constexpr static integer index[NDIM][NDIM] = { { 0, 1, 2 }, { 1, 3, 4 }, {
			2, 4, 5 } };
	simd_vector data[SIZE];
public:
	tensor() = default;
	tensor(const tensor&) = default;
	tensor(tensor&&) = default;
	~tensor() = default;
	tensor& operator=(const tensor&) = default;
	tensor& operator=(tensor &&) = default;
	const simd_vector& operator()(integer i, integer j) const {
		return data[index[i][j]];
	}
	simd_vector& operator()(integer i, integer j) {
		return data[index[i][j]];
	}
	tensor& operator*=(const scalar& other) {
		for (integer i = 0; i != NDIM; ++i) {
			data[i] *= other;
		}
		return *this;
	}
	tensor& operator/=(const scalar& other) {
		const auto other_inv = one / other;
		return operator*=(other_inv);
	}
	tensor operator*(const scalar& other) const {
		tensor c(*this);
		for (integer i = 0; i != NDIM; ++i) {
			c.data[i] *= other;
		}
		return c;
	}
	tensor operator/(const scalar& other) const {
		const auto other_inv = one / other;
		return operator*(other_inv);
	}

	tensor operator*(const tensor& other) const {
		tensor c;
		for (integer i = 0; i != NDIM; ++i) {
			for (integer j = 0; j <= i; ++j) {
				simd_vector& sum = c(i, j);
				sum = +(*this)(i, 0) * other(0, j);
				sum += (*this)(i, 1) * other(1, j);
				sum += (*this)(i, 2) * other(2, j);
			}
		}
		return c;
	}

	tensor& operator*=(const tensor& other) {
		auto c = (*this) * other;
		(*this) = std::move(c);
		return *this;
	}

	vector operator*(const vector& other) const {
		vector c;
		for (integer i = 0; i != NDIM; ++i) {
			c(i) = +(*this)(i, 0) * other(0);
			c(i) += (*this)(i, 1) * other(1);
			c(i) += (*this)(i, 2) * other(1);
		}
		return c;
	}
};

constexpr integer tensor::SIZE;
constexpr integer tensor::index[NDIM][NDIM];

union z4_t {
	static constexpr integer NFIELD = 38;
	typedef struct {
		scalar alpha;
		tensor gamma;
		scalar At1;
		scalar At2;
		scalar Ep;
		scalar Em;
		scalar Gp;
		scalar Gm;
		vector S;
		tensor Dt1;
		tensor Dt2;
		tensor Lp;
		tensor Lm;
	} z4_char;
	struct {
		scalar alpha;
		tensor gamma;
		scalar Theta;
		vector A;
		vector Z;
		tensor K;
		tensor D[3];
	} f;
	simd_vector array[NFIELD];
	z4_t(const z4_t& other) {
		for (integer f = 0; f != NFIELD; ++f) {
			array[f] = other.array[f];
		}
	}
	static void suppress_char_dir(z4_char& c, integer pos) {
		if (pos > 0) {
			c.Ep = c.Gp = zero;
			for (integer i = 0; i != NDIM; ++i) {
				for (integer j = i; j != NDIM; ++j) {
					c.Lp(i, j) = zero;
				}
			}
		} else {
			c.Em = c.Gm = zero;
			for (integer i = 0; i != NDIM; ++i) {
				for (integer j = i; j != NDIM; ++j) {
					c.Lm(i, j) = zero;
				}
			}

		}
	}
	void to_char_inv(const z4_char& c, integer dim) {
		integer ti1, ti2;
		vector n;
		if (dim == 0) {
			ti1 = 1;
			ti2 = 2;
			n(0) = one;
			n(1) = n(2) = zero;
		} else if (dim == 1) {
			ti1 = 0;
			ti2 = 2;
			n(1) = one;
			n(0) = n(2) = zero;
		} else {
			ti1 = 0;
			ti2 = 1;
			n(2) = one;
			n(1) = n(0) = zero;
		}
		auto& u = f;
		u.gamma = f.gamma;
		u.alpha = f.alpha;

		u.A(dim) = half * (c.Gp - c.Gm);
		u.A(ti1) = c.At1;
		u.A(ti2) = c.At2;
		u.Theta = half * (c.Em + c.Ep);
		const auto Vn = half * (c.Ep - c.Em);
		const auto K = half * (c.Gp + c.Gm) + two * u.Theta;
		for (integer i = 0; i != NDIM; ++i) {
			for (integer j = i; j != NDIM; ++j) {
				u.D[ti1](i, j) = c.Dt1(i, j);
				u.D[ti2](i, j) = c.Dt2(i, j);
				u.K(i, j) = half * (c.Lp(i, j) + c.Lm(i, j));
				u.K(i, j) += n(i) * n(j) * K;
				u.D[dim](i, j) = half * (c.Lp(i, j) - c.Lm(i, j));
				u.D[dim](i, j) -= half * c.S(i) * delta(dim, j);
				u.D[dim](i, j) -= half * c.S(j) * delta(dim, i);
			}
		}
		u.D[dim](dim, dim) = zero;
		scalar Dstar = zero;
		for (integer i = 0; i != NDIM; ++i) {
			Dstar += u.D[dim](i, i);
		}
		u.D[dim](dim, dim) = u.A(dim) - Dstar + two * Vn - c.S(dim);
		vector D;
		vector E;
		for (integer k = 0; k != NDIM; ++k) {
			D(k) = E(k) = zero;
			for (integer i = 0; i != NDIM; ++i) {
				D(k) += u.D[k](i, i);
				E(k) += u.D[i](k, i);
			}
		}
		for (integer i = 0; i != NDIM; ++i) {
			u.Z(i) = half * (u.A(i) + D(i) - two * E(i) - c.S(i));
		}
	}
	z4_char to_char(integer dim) {
		z4_char c;
		integer ti1, ti2;
		vector n;
		if (dim == 0) {
			ti1 = 1;
			ti2 = 2;
			n(0) = one;
			n(1) = n(2) = zero;
		} else if (dim == 1) {
			ti1 = 0;
			ti2 = 2;
			n(1) = one;
			n(0) = n(2) = zero;
		} else {
			ti1 = 0;
			ti2 = 1;
			n(2) = one;
			n(1) = n(0) = zero;
		}
		scalar K;
		vector E;
		vector D;
		vector V;
		tensor lambda_k;
		compute_traces(K, E, D);
		auto lambda_n = zero;
		for (integer i = 0; i != NDIM; ++i) {
			V(i) = D(i) - E(i) - f.Z(i);
		}
		for (integer i = 0; i != NDIM; ++i) {
			for (integer j = i; j != NDIM; ++j) {
				auto& lam = lambda_k(i, j);
				lam = 0.0;
				lam += f.D[dim](i, j);
				if (i == dim) {
					lam += half * (f.A(j) - D(j));
					lam += V(j);
				}
				if (j == dim) {
					lam += half * (f.A(i) - D(i));
					lam += V(i);
				}
			}
		}
		for (integer i = 0; i != NDIM; ++i) {
			lambda_n += lambda_k(i, i);
		}

		scalar Vn = V(dim);
		z4_t U = *this;
		const auto& u = U.f;
		c.alpha = u.alpha;
		c.gamma = u.gamma;
		c.At1 = u.A(ti1);
		c.At2 = u.A(ti2);
		c.Ep = u.Theta + Vn;
		c.Em = u.Theta - Vn;
		const auto term1 = K - two * u.Theta;
		simd_vector term2 = u.A(dim);
		for (integer k = 0; k != NDIM; ++k) {
			c.S(k) = u.A(k) - D(k) + two * V(k);
		}
		c.Gp = term1 + term2;
		c.Gm = term1 - term2;
		for (integer i = 0; i != NDIM; ++i) {
			for (integer j = i; j != NDIM; ++j) {
				c.Dt1(i, j) = u.D[ti1](i, j);
				c.Dt2(i, j) = u.D[ti2](i, j);
				const auto term1 = u.K(i, j) - n(i) * n(j) * K;
				const auto term2 = lambda_k(i, j) - n(i) * n(j) * lambda_n;
				c.Lp(i, j) = term1 + term2;
				c.Lm(i, j) = term1 - term2;
			}
		}
		return c;
	}
	simd_vector compute_sqrt_gamma() const {
		return one + half * (f.gamma(0, 0) + f.gamma(1, 1) + f.gamma(2, 2));
	}
	void compute_traces(scalar& K, vector& E, vector& D) const {
		K = zero;
		for (integer i = 0; i != NDIM; ++i) {
			auto& e = E(i);
			auto& d = D(i);
			e = zero;
			d = zero;
			K += f.K(i, i);
			for (integer j = 0; j != NDIM; ++j) {
				e += f.D[j](i, j);
				d += f.D[i](j, j);
			}
		}
	}

	void flux(z4_t& F, simd_vector& a, integer k,
			const std::array<simd_vector, NDIM>& X) const {
		scalar K;
		vector E;
		vector D;
		vector V;
		tensor lambda_k;
		compute_traces(K, E, D);
		for (integer i = 0; i != NDIM; ++i) {
			V(i) = D(i) - E(i) - f.Z(i);
		}
		for (integer i = 0; i != NDIM; ++i) {
			for (integer j = i; j != NDIM; ++j) {
				auto& lam = lambda_k(i, j);
				lam = 0.0;
				lam += f.D[k](i, j);
				if (i == k) {
					lam += half * (f.A(j) - D(j));
					lam += V(j);
				}
				if (j == k) {
					lam += half * (f.A(i) - D(i));
					lam += V(i);
				}
			}
		}

		std::fill(F.array, F.array + NFIELD, zero);

		const auto Q = K - two * f.Theta;
		F.f.Z(k) += (K - f.Theta);

		F.f.Theta += (D(k) - E(k) - f.Z(k));

		F.f.A(k) += (Q - eta * f.alpha);

		for (integer i = 0; i != NDIM; ++i) {
			F.f.Z(i) -= f.K(i, k);
		}

		for (integer i = 0; i != NDIM; ++i) {
			for (integer j = i; j != NDIM; ++j) {
				F.f.K(i, j) += lambda_k(i, j);
			}
		}

		for (integer i = 0; i != NDIM; ++i) {
			for (integer j = i; j != NDIM; ++j) {
				F.f.D[k](i, j) += f.K(i, j);
				F.f.D[k](i, j) -= eta * half * f.gamma(i, j);
			}
		}
		a = one;

	}
	void source(z4_t& S, const scalar& tau, const vector& Si, const tensor& Sij,
			const std::array<simd_vector, NDIM>& X) const {
		std::fill(S.array, S.array + NFIELD, simd_vector(0));
		static const simd_vector eight_pi(8.0 * M_PI);
		static const simd_vector four_pi(4.0 * M_PI);
		scalar TrS = zero;
		scalar TrK = zero;
		vector E;
		vector D;

		compute_traces(TrK, E, D);

		for (integer i = 0; i != NDIM; ++i) {
			TrS += Sij(i, i);
		}

		const auto Q = (TrK - two * f.Theta);

		S.f.alpha = -Q;

		for (integer i = 0; i != NDIM; ++i) {
			S.f.A(i) = -eta * f.A(i);
		}

		for (integer i = 0; i != NDIM; ++i) {
			for (integer j = i; j != NDIM; ++j) {
				S.f.gamma(i, j) = -two * f.K(i, j);
			}
		}

		for (integer i = 0; i != NDIM; ++i) {
			for (integer j = i; j != NDIM; ++j) {
				for (integer k = 0; k != NDIM; ++k) {
					S.f.D[k](i, j) = -eta * f.D[k](i, j);
				}
			}
		}

		for (integer i = 0; i != NDIM; ++i) {
			S.f.K(i, i) -= kap1 * (one + kap2) * f.Theta;
		}

		for (integer i = 0; i != NDIM; ++i) {
			S.f.Z(i) -= kap1 * f.Z(i);
		}
		S.f.Theta -= kap1 * (one + half * kap2) * f.Theta;

		for (integer i = 0; i != NDIM; ++i) {
			for (integer j = i; j != NDIM; ++j) {
				auto tmp = zero;
				tmp -= eight_pi * Sij(i, j);
				tmp += four_pi * (TrS - tau) * delta(i, j);
				S.f.K(i, j) += tmp;
			}
		}

		for (integer i = 0; i != NDIM; ++i) {
			S.f.Z(i) -= eight_pi * Si(i);
		}

		S.f.Theta -= eight_pi * tau;

	}

	z4_t() :
			f() {
	}
}
;

}

gr::z4_t& grhd::get_z4(vector_type& ref) {
	return *reinterpret_cast<gr::z4_t*>(ref.data() + srhd::NF);
}

const gr::z4_t& grhd::get_z4(const vector_type& ref) {
	return *reinterpret_cast<const gr::z4_t*>(ref.data() + srhd::NF);
}

void grhd::initial_value(vector_type& u, const std::array<simd_vector, NDIM>& x,
		real) {
	auto& z4 = get_z4(u);
	for (integer i = 0; i != simd_len; ++i) {
		real r = x[XDIM][i] * x[XDIM][i];
		r += x[YDIM][i] * x[YDIM][i];
		r += x[ZDIM][i] * x[ZDIM][i];
		r = sqrt(r);
		const real r0 = 1.0;
		const real v0 = r0 * r0 * r0 * 4.0 / 3.0 * M_PI;
		const real d0 = 1.0e-3 / v0;
		u[D_i][i] = r < r0 ? d0 : 1.0e-6 * d0;
		u[tau_i][i] = u[D_i][i];
		u[Sx_i][i] = u[Sy_i][i] = u[Sz_i][i] = 0.0;
		for (integer j = 0; j != gr::z4_t::NFIELD; ++j) {
			z4.array[j][i] = 0.0;
		}
		z4.f.alpha[i] = 0.;
		z4.f.gamma(0, 0)[i] = 0.;
		z4.f.gamma(1, 1)[i] = 0.;
		z4.f.gamma(2, 2)[i] = 0.;
	}
}

void grhd::set_outflow(vector_type& u, const geo::face& face, real dx) {
	const auto d = face.get_dimension();
	const simd_vector s(2 * face.get_side() - 1);
	if (face.get_side() > 0) {
		u[Sx_i + d] = max(u[Sx_i + d], zero);
	} else {
		u[Sx_i + d] = min(u[Sx_i + d], zero);
	}
	auto& u_z4 = get_z4(u);
	auto v = u_z4;
	auto c = u_z4.to_char(d);
	u_z4.suppress_char_dir(c, -2 * face.get_side() + 1);
	u_z4.to_char_inv(c, d);
	u_z4.f.alpha = 0.0;
	for (integer i = 0; i != NDIM; ++i) {
		for (integer j = i; j != NDIM; ++j) {
			u_z4.f.gamma(i, j) = 0.0;
		}
	}
}

void grhd::physical_flux(vector_type& f, simd_vector& a, const vector_type& u,
		const vector_type& v, integer dim,
		const std::array<simd_vector, NDIM>& X, real t) {

	const auto& v_z4 = get_z4(v);

	const auto& alpha = v_z4.f.alpha;

	const auto& D = u[D_i];
	const auto& tau = u[tau_i];
	const auto& Sx = u[Sx_i];
	const auto& Sy = u[Sy_i];
	const auto& Sz = u[Sz_i];

	const auto& rho = v[rho_i];
	const auto& eps = v[eps_i];
	const auto& uel = v[ux_i + dim];
	const auto& ux = v[ux_i];
	const auto& uy = v[uy_i];
	const auto& uz = v[uz_i];

	const auto u2 = ux * ux + uy * uy + uz * uz;
	const auto W = sqrt(one + u2);
	const auto Winv = one / W;
	const auto vel = uel * Winv;
	const auto v2 = u2 * Winv * Winv;
	const auto p = (fgamma - one) * rho * eps;
	const auto h = rho * (one + eps) + p;
	const auto hinv = one / h;
	const auto c2 = fgamma * p * hinv;
	const auto vel2 = vel * vel;
	const auto ss = vel;
	if (t > 4.0 * opts.xscale) {
		f[D_i] = D * ss;
		f[Sx_i] = Sx * ss;
		f[Sy_i] = Sy * ss;
		f[Sz_i] = Sz * ss;
		f[Sx_i + dim] += p;
		f[tau_i] = (tau + p) * ss;
	}
	const auto onemv2c2 = one - v2 * c2;
	const auto onemv2c2inv = one / onemv2c2;
	const auto onemc2 = one - c2;
	const auto a1 = onemv2c2inv * onemc2 * vel;
	const auto a2 = onemv2c2inv * sqrt(c2 * Winv * (onemv2c2 - onemc2 * vel2));
	a = (one + alpha) * (abs(a1) + a2);

	if (t > 4.0 * opts.xscale) {
		const auto speed = one + v_z4.f.alpha;
		for (integer field = 0; field != srhd::NF; ++field) {
			f[field] *= speed;
		}
	}
}

void grhd::to_con(vector_type& u, const vector_type& v) {

	const auto& v_z4 = get_z4(v);

	const auto sqrt_gamma = v_z4.compute_sqrt_gamma();

	auto& D = u[D_i];
	auto& tau = u[tau_i];
	auto& Sx = u[Sx_i];
	auto& Sy = u[Sy_i];
	auto& Sz = u[Sz_i];

	const auto& rho = v[rho_i];
	const auto& eps = v[eps_i];
	const auto& ux = v[ux_i];
	const auto& uy = v[uy_i];
	const auto& uz = v[uz_i];

	const auto u2 = ux * ux + uy * uy + uz * uz;
	const auto W = sqrt(one + u2);
	const auto p = (fgamma - 1.0) * rho * eps;
	const auto h = rho * (one + fgamma * eps);
	const auto hW = W * h;
	D = W * rho * sqrt_gamma;
	tau = (hW * W - p) * sqrt_gamma;
	Sx = hW * ux * sqrt_gamma;
	Sy = hW * uy * sqrt_gamma;
	Sz = hW * uz * sqrt_gamma;

	for (integer f = srhd::NF; f != NF; ++f) {
		u[f] = v[f];
	}
}

void grhd::to_prim(vector_type& v, const vector_type& u) {
	static const simd_vector k((fgamma - one) / fgamma);

	const auto& u_z4 = get_z4(u);

	const auto sqrt_gamma = u_z4.compute_sqrt_gamma();

	const auto D = u[D_i] / sqrt_gamma;
	auto tau = u[tau_i] / sqrt_gamma;
	const auto Sx = u[Sx_i] / sqrt_gamma;
	const auto Sy = u[Sy_i] / sqrt_gamma;
	const auto Sz = u[Sz_i] / sqrt_gamma;

	auto& rho = v[rho_i];
	auto& eps = v[eps_i];
	auto& ux = v[ux_i];
	auto& uy = v[uy_i];
	auto& uz = v[uz_i];

	simd_vector W, vel, f, dfdv;
	simd_vector S(sqrt(Sx * Sx + Sy * Sy + Sz * Sz));

	tau = max(tau, sqrt(S * S + D * D));

	vel = S / D;
	const auto Sinv = S.inv_or_zero();
	vel = min(abs(vel), simd_vector(0.999));
	for (int i = 0; i < 8; ++i) {
		W = sqrt(one / (one - vel * vel));
		const simd_vector Winv(one / W);
		const simd_vector W2inv(Winv * Winv);
		f = (one - k * W2inv) * S + vel * k * D * Winv - vel * tau;
		dfdv = k * (two * S * vel + D * (one - two * vel * vel) * W) - tau;
		const auto v0 = vel;
		vel = vel - f / dfdv;
		const simd_vector b(0.1);
		vel = max(b * v0, vel);
		vel = min(one - b * (one - v0), vel);
	}
	W = sqrt(one / (one - vel * vel));
	const auto SinvvelW = Sinv * vel * W;
	rho = D / W;
	eps = (tau / rho - W * W) / (fgamma * W * W - (fgamma - one));
	ux = Sx * SinvvelW;
	uy = Sy * SinvvelW;
	uz = Sz * SinvvelW;
	eps = max(eps, zero);

	for (integer f = srhd::NF; f != NF; ++f) {
		v[f] = u[f];
	}
}

void grhd::to_output(vector_type& v, const vector_type& u) {
	to_prim(v, u);
	auto& ux = v[ux_i];
	auto& uy = v[uy_i];
	auto& uz = v[uz_i];
	const auto u2 = ux * ux + uy * uy + uz * uz;
	const auto W = sqrt(one + u2);
	const auto Winv = one / W;
	ux *= Winv;
	uy *= Winv;
	uz *= Winv;
	for (integer f = srhd::NF; f != NF; ++f) {
		v[f] = u[f];
	}
}

void grhd::to_fluxes(vector_type& flux, simd_vector& a, vector_type& vl,
		vector_type& vr, integer dim, const std::array<simd_vector, NDIM>& X,
		real t) {

	vector_type ur, ul, fr, fl;
	simd_vector afr, afl, agl, agr;

	to_con(ur, vr);
	to_con(ul, vl);

	auto& fl_z4 = get_z4(fl);
	auto& fr_z4 = get_z4(fr);
	auto& ul_z4 = get_z4(ul);
	auto& ur_z4 = get_z4(ur);
	physical_flux(fr, afr, ur, vr, dim, X, t);
	physical_flux(fl, afl, ul, vl, dim, X, t);
	ur_z4.flux(fr_z4, agr, dim, X);
	ul_z4.flux(fl_z4, agl, dim, X);
	const auto af = max(afr, afl);
	const auto ag = max(agr, agl);
	if (t > 4.0 * opts.xscale) {
		for (integer f = 0; f != srhd::NF; ++f) {
			flux[f] = half * (fr[f] + fl[f]);
			flux[f] -= half * af * (ur[f] - ul[f]);
		}
	}
	for (integer f = srhd::NF; f != NF; ++f) {
		const integer n = f - srhd::NF;
		flux[f] = half * (fr[f] + fl[f]);
		flux[f] -= half * ag * (ur[f] - ul[f]);
	}
	a = max(ag, af);
}

bool grhd::refinement_test(integer level,
		const std::array<simd_vector, NDIM>& x, const vector_type& u,
		const std::array<vector_type, NDIM>& dudx) {
	const auto max_level = opts.max_level;
	return level < max_level;
}

void grhd::explicit_source(vector_type& s, const vector_type& u,
		const vector_type& v, const std::array<simd_vector, NDIM>& X, real t) {
	std::fill(s.begin(), s.end(), simd_vector(0));

	gr::scalar tau0;
	gr::vector Si;
	gr::tensor Sij;

	const auto& tau = u[tau_i];
	const auto& Sx = u[Sx_i];
	const auto& Sy = u[Sy_i];
	const auto& Sz = u[Sz_i];

	const auto& ux = v[ux_i];
	const auto& uy = v[uy_i];
	const auto& uz = v[uz_i];

	const auto u2 = ux * ux + uy * uy + uz * uz;
	const auto W = sqrt(one + u2);
	const auto Winv = one / W;

	tau0 = tau;
	Si(XDIM) = Sx;
	Si(YDIM) = Sy;
	Si(ZDIM) = Sz;
	Sij(XDIM, XDIM) = Sx * ux * Winv;
	Sij(XDIM, YDIM) = Sx * uy * Winv;
	Sij(XDIM, ZDIM) = Sx * uz * Winv;
	Sij(YDIM, YDIM) = Sy * uy * Winv;
	Sij(YDIM, ZDIM) = Sy * uz * Winv;
	Sij(ZDIM, ZDIM) = Sz * uz * Winv;

	auto& s_z4 = get_z4(s);
	auto& u_z4 = get_z4(u);

	u_z4.source(s_z4, tau0, Si, Sij, X);
	if (t > 4.0 * opts.xscale) {
		for (integer i = 0; i != NDIM; ++i) {
			s[tau_i] -= Si(i) * u_z4.f.A(i);
			s[Sx_i + i] -= u_z4.f.A(i) * tau0;
			for (integer j = 0; j != NDIM; ++j) {
				s[tau_i] += Sij(i, j) * u_z4.f.K(i, j);
				for (integer k = 0; k != NDIM; ++k) {
					s[Sx_i + k] += u_z4.f.D[k](i, j) * Sij(i, j);
				}
			}
		}
	}

}

void grhd::implicit_source(vector_type& s, const vector_type& u, real) {
	std::fill(s.begin(), s.end(), simd_vector(0));

}

std::vector<std::string> grhd::field_names() {
	std::vector<std::string> names(NF);
	names[D_i] = "rho";
	names[tau_i] = "eps";
	names[Sx_i] = "vx";
	names[Sy_i] = "vy";
	names[Sz_i] = "vz";
	names[srhd::NF + 0] = "alpha";
	names[srhd::NF + 1] = "gamma_xx";
	names[srhd::NF + 2] = "gamma_xy";
	names[srhd::NF + 3] = "gamma_xz";
	names[srhd::NF + 4] = "gamma_yy";
	names[srhd::NF + 5] = "gamma_yz";
	names[srhd::NF + 6] = "gamma_zz";
	names[srhd::NF + 7] = "Theta";
	names[srhd::NF + 8] = "A_x";
	names[srhd::NF + 9] = "A_y";
	names[srhd::NF + 10] = "A_z";
	names[srhd::NF + 11] = "Z_x";
	names[srhd::NF + 12] = "Z_y";
	names[srhd::NF + 13] = "Z_z";
	names[srhd::NF + 14] = "K_xx";
	names[srhd::NF + 15] = "K_xy";
	names[srhd::NF + 16] = "K_xz";
	names[srhd::NF + 17] = "K_yy";
	names[srhd::NF + 18] = "K_yz";
	names[srhd::NF + 19] = "K_zz";
	names[srhd::NF + 20] = "Dx_xx";
	names[srhd::NF + 21] = "Dx_xy";
	names[srhd::NF + 22] = "Dx_xz";
	names[srhd::NF + 23] = "Dx_yy";
	names[srhd::NF + 24] = "Dx_yz";
	names[srhd::NF + 25] = "Dx_zz";
	names[srhd::NF + 26] = "Dy_xx";
	names[srhd::NF + 27] = "Dy_xy";
	names[srhd::NF + 28] = "Dy_xz";
	names[srhd::NF + 29] = "Dy_yy";
	names[srhd::NF + 30] = "Dy_yz";
	names[srhd::NF + 31] = "Dy_zz";
	names[srhd::NF + 32] = "Dz_xx";
	names[srhd::NF + 33] = "Dz_xy";
	names[srhd::NF + 34] = "Dz_xz";
	names[srhd::NF + 35] = "Dz_yy";
	names[srhd::NF + 36] = "Dz_yz";
	names[srhd::NF + 37] = "Dz_zz";
	return names;
}

