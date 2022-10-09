// Estimating temporally variable selection intensity from ancient DNA data
// Zhangyi He, Xiaoyang Dai, Wenyang Lyu, Mark Beaumont, Feng Yu

// version 1.5 (PMMH)
// Phenotypes controlled by a single gene
// Non-constant natural selection and non-constant demographic histories

// Input: genotype likelihoods
// Output: posteriors for the selection coefficient, dominance parameter and the genotype frequency trajectories of the population

// C functions

#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]
#include <math.h>

using namespace Rcpp;
using namespace std;
// using namespace arma;

/********** WFM **********/
// Calculate the fitness matrix for the single-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::dmat calculateFitnessMat_arma(const double& sel_cof, const double& dom_par) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // calculate the fitness matrix
  arma::dmat fts_mat = arma::ones<arma::dmat>(2, 2);
  // fts_mat(0, 0) = 1;
  fts_mat(1, 0) = 1 + sel_cof * dom_par;
  fts_mat(0, 1) = 1 + sel_cof * dom_par;
  fts_mat(1, 1) = 1 + sel_cof;

  // return the fitness matrix for the Wright-Fisher model
  return fts_mat;
}

// Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher model with selection
// [[Rcpp::export]]
List simulateWFM_arma(const arma::dmat& fts_mat, const arma::icolvec& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the allele and genotype frequency trajectories
  arma::drowvec mut_frq_pth = arma::zeros<arma::drowvec>(arma::uword(lst_gen - int_gen) + 1);
  arma::dmat gen_frq_pth = arma::zeros<arma::dmat>(3, arma::uword(lst_gen - int_gen) + 1);

  // initialise the allele and genotype frequencies in generation 0
  arma::dcolvec ale_frq = {1 - int_frq, int_frq};
  mut_frq_pth(0) = ale_frq(1);
  arma::dmat gen_frq = fts_mat % (ale_frq * ale_frq.t()) / arma::as_scalar(ale_frq.t() * fts_mat * ale_frq);
  gen_frq = arma::diagmat(gen_frq) + 2 * arma::trimatu(gen_frq, 1);
  gen_frq_pth.col(0) = gen_frq(arma::trimatu_ind(arma::size(gen_frq)));

  // simulate the allele and genotype frequency trajectories
  for (arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    // calculate the sampling probabilities
    arma::dcolvec prob = ale_frq;
    prob = ale_frq % (fts_mat * ale_frq) / arma::as_scalar(ale_frq.t() * fts_mat * ale_frq);

    // proceed the Wright-Fisher sampling
    ale_frq(0) = R::rbinom(2 * pop_siz(k), prob(0)) / 2 / pop_siz(k);
    ale_frq(1) = 1 - ale_frq(0);
    mut_frq_pth(k) = ale_frq(1);

    gen_frq = fts_mat % (ale_frq * ale_frq.t()) / arma::as_scalar(ale_frq.t() * fts_mat * ale_frq);
    gen_frq = arma::diagmat(gen_frq) + 2 * arma::trimatu(gen_frq, 1);
    gen_frq_pth.col(k) = gen_frq(arma::trimatu_ind(arma::size(gen_frq)));
  }

  // return the mutant allele and genotype frequency trajectories under the Wright-Fisher model
  return List::create(Named("mut_frq_pth", mut_frq_pth),
                      Named("gen_frq_pth", gen_frq_pth));
}
/*************************/


/********** WFD **********/
// Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::drowvec simulateWFD_arma(const double& sel_cof, const double& dom_par, const arma::icolvec& pop_siz, const int& ref_siz, const double& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // rescale the selection coefficient
  double scl_sel_cof = 2 * ref_siz * sel_cof;

  // calculate the ratio of the population size to the reference population size
  arma::dcolvec siz_rto = arma::zeros<arma::dcolvec>(arma::uword(lst_gen - int_gen) * ptn_num);
  for (arma::uword k = 0; k < arma::uword(lst_gen - int_gen); k++) {
    arma::dcolvec siz_rto_tmp = arma::zeros<arma::dcolvec>(ptn_num);
    siz_rto_tmp.fill(pop_siz(k) / double(ref_siz));
    siz_rto.subvec(k * ptn_num, (k + 1) * ptn_num - 1) = siz_rto_tmp;
  }

  // declare delta t
  double dt = 1.0 / (2 * ref_siz) / ptn_num;
  // declare delta W
  arma::drowvec dW = pow(dt, 0.5) * arma::randn<arma::drowvec>(arma::uword(lst_gen - int_gen) * ptn_num);

  // declare the mutant allele frequency trajectory
  arma::drowvec mut_frq_pth = arma::zeros<arma::drowvec>(arma::uword(lst_gen - int_gen) * ptn_num + 1);

  // initialise the mutant allele frequency in generation 0
  mut_frq_pth(0) = int_frq;

  // simulate the mutant allele frequency trajectory
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient
    double mu = scl_sel_cof * mut_frq_pth(t - 1) * (1 - mut_frq_pth(t - 1)) * (dom_par + (1 - 2 * dom_par) * mut_frq_pth(t - 1));

    // calculate the diffusion coefficient
    double sigma = pow(mut_frq_pth(t - 1) * (1 - mut_frq_pth(t - 1)) / siz_rto(t - 1), 0.5);

    // proceed the Euler-Maruyama scheme
    mut_frq_pth(t) = mut_frq_pth(t - 1) + mu * dt + sigma * dW(t - 1);

    // remove the noise from the numerical techniques
    if (mut_frq_pth(t) < 0) {
      mut_frq_pth(t) = 0;
    }
    if (mut_frq_pth(t) > 1) {
      mut_frq_pth(t) = 1;
    }
  }

  // return the mutant allele frequency trajectory under the Wright-Fisher diffusion
  return mut_frq_pth;
}
/*************************/


/********** BPF **********/
// Calculate the genotype frequencies in the adult with the mutant allele frequency
// [[Rcpp::export]]
arma::dmat calculateGenoFrq_arma(const double& sel_cof, const double& dom_par, const arma::dcolvec& mut_frq) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat gen_frq = arma::zeros<arma::dmat>(mut_frq.n_elem, 3);
  gen_frq.col(0) = pow(1 - mut_frq, 2);
  gen_frq.col(1) = (1 + dom_par * sel_cof) * 2 * (mut_frq % (1 - mut_frq));
  gen_frq.col(2) = (1 + sel_cof) * pow(mut_frq, 2);
  gen_frq = arma::normalise(gen_frq, 1, 1);

  return gen_frq.t();
}

// Calculate the particle weights
// [[Rcpp::export]]
arma::dcolvec calculateParticleWeight_arma(const arma::dmat& gen_lik, const double& sel_cof, const double& dom_par, const arma::dcolvec& mut_frq) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dmat gen_frq = arma::zeros<arma::dmat>(mut_frq.n_elem, 3);
  gen_frq.col(0) = pow(1 - mut_frq, 2);
  gen_frq.col(1) = (1 + dom_par * sel_cof) * 2 * (mut_frq % (1 - mut_frq));
  gen_frq.col(2) = (1 + sel_cof) * pow(mut_frq, 2);
  gen_frq = arma::normalise(gen_frq, 1, 1);

  arma::dcolvec prob = arma::prod(gen_frq * gen_lik, 1);

  return prob;
}

// Run the bootstrap particle filter
// [[Rcpp::export]]
List runBPF_arma(const arma::dcolvec& sel_cof, const double& dom_par, const arma::icolvec& pop_siz, const int& ref_siz, const int& evt_gen, const arma::dmat& raw_smp, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::irowvec smp_gen = arma::unique(arma::conv_to<arma::irowvec>::from(raw_smp.row(0)));
  arma::field<arma::dmat> gen_lik(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    arma::ucolvec idx = arma::find(raw_smp.row(0) == smp_gen(k));
    gen_lik(k) = raw_smp.submat(1, idx.min(), 3, idx.max());
  }

  double log_lik = 0;
  double lik = exp(log_lik);

  arma::uword evt_idx = arma::as_scalar(arma::max(arma::find(smp_gen <= evt_gen))) + 1;

  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dmat mut_frq_pth = arma::zeros<arma::dmat>(pcl_num, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1);
  arma::dmat mut_frq_pre = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dmat mut_frq_pst = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dcube gen_frq_pre = arma::zeros<arma::dcube>(3, pcl_num, smp_gen.n_elem);
  arma::dcube gen_frq_pst = arma::zeros<arma::dcube>(3, pcl_num, smp_gen.n_elem);

  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dcolvec mut_frq_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat gen_frq_tmp = arma::zeros<arma::dmat>(3, pcl_num);

  // before the event of interest
  // initialise the particles
  cout << "generation: " << smp_gen(0) << endl;
  mut_frq_tmp = arma::randu<arma::dcolvec>(pcl_num);
  gen_frq_tmp = calculateGenoFrq_arma(sel_cof(0), dom_par, mut_frq_tmp);
  wght_tmp = calculateParticleWeight_arma(gen_lik(0), sel_cof(0), dom_par, mut_frq_tmp);

  if (arma::sum(wght_tmp) > 0) {
    // arma::dcolvec prob = arma::normalise(wght_tmp, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, wght_tmp);

    lik = lik * arma::mean(wght_tmp);
    wght.col(0) = wght_tmp;
    mut_frq_pre.col(0) = mut_frq_tmp;
    mut_frq_pst.col(0) = mut_frq_tmp.elem(indx);
    gen_frq_pre.slice(0) = gen_frq_tmp;
    gen_frq_pst.slice(0) = gen_frq_tmp.cols(indx);
  } else {
    lik = 0;
    wght.shed_cols(0, smp_gen.n_elem - 1);
    mut_frq_pre.shed_cols(0, smp_gen.n_elem - 1);
    mut_frq_pst.shed_cols(0, smp_gen.n_elem - 1);
    gen_frq_pre.shed_slices(0, smp_gen.n_elem - 1);
    gen_frq_pst.shed_slices(0, smp_gen.n_elem - 1);

    return List::create(Named("lik", lik),
                        Named("wght", wght),
                        Named("mut_frq_pth", mut_frq_pth),
                        Named("mut_frq_pre_resmp", mut_frq_pre),
                        Named("mut_frq_pst_resmp", mut_frq_pst),
                        Named("gen_frq_pre_resmp", gen_frq_pre),
                        Named("gen_frq_pst_resmp", gen_frq_pst));
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < evt_idx; k++) {
    cout << "generation: " << smp_gen(k) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    mut_frq_tmp = mut_frq_pst.col(k - 1);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::drowvec path = simulateWFD_arma(sel_cof(0), dom_par, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, mut_frq_tmp(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      mut_frq_pth.submat(i, (smp_gen(k - 1) - smp_gen(0)) * ptn_num, i, (smp_gen(k) - smp_gen(0)) * ptn_num) = path;
      mut_frq_tmp(i) = arma::as_scalar(path.tail(1));
    }
    gen_frq_tmp = calculateGenoFrq_arma(sel_cof(0), dom_par, mut_frq_tmp);
    wght_tmp = calculateParticleWeight_arma(gen_lik(k), sel_cof(0), dom_par, mut_frq_tmp);

    if (arma::sum(wght_tmp) > 0) {
      // arma::dcolvec prob = arma::normalise(wght_tmp, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, wght_tmp);

      lik = lik * arma::mean(wght_tmp);
      wght.col(k) = wght_tmp;
      mut_frq_pre.col(k) = mut_frq_tmp;
      mut_frq_pst.col(k) = mut_frq_tmp.elem(indx);
      mut_frq_pth = mut_frq_pth.rows(indx);
      gen_frq_pre.slice(k) = gen_frq_tmp;
      gen_frq_pst.slice(k) = gen_frq_tmp.cols(indx);
    } else {
      lik = 0;
      wght.shed_cols(k, smp_gen.n_elem - 1);
      mut_frq_pre.shed_cols(k, smp_gen.n_elem - 1);
      mut_frq_pst.shed_cols(k, smp_gen.n_elem - 1);
      gen_frq_pre.shed_slices(k, smp_gen.n_elem - 1);
      gen_frq_pst.shed_slices(k, smp_gen.n_elem - 1);

      return List::create(Named("lik", lik),
                          Named("wght", wght),
                          Named("mut_frq_pth", mut_frq_pth),
                          Named("mut_frq_pre_resmp", mut_frq_pre),
                          Named("mut_frq_pst_resmp", mut_frq_pst),
                          Named("gen_frq_pre_resmp", gen_frq_pre),
                          Named("gen_frq_pst_resmp", gen_frq_pst));
    }
  }

  // after the event of interest
  // run the bootstrap particle filter
  if (evt_idx < smp_gen.n_elem) {
    if (smp_gen(evt_idx - 1) != evt_gen) {
      cout << "generation: " << evt_gen << endl;
      mut_frq_tmp = mut_frq_pst.col(evt_idx - 1);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateWFD_arma(sel_cof(0), dom_par, pop_siz.subvec(smp_gen(evt_idx - 1), evt_gen), ref_siz, mut_frq_tmp(i), smp_gen(evt_idx - 1), evt_gen, ptn_num);
        mut_frq_pth.submat(i, (smp_gen(evt_idx - 1) - smp_gen(0)) * ptn_num, i, (evt_gen - smp_gen(0)) * ptn_num) = path;
        mut_frq_tmp(i) = arma::as_scalar(path.tail(1));
      }

      cout << "generation: " << smp_gen(evt_idx) << endl;
      wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateWFD_arma(sel_cof(1), dom_par, pop_siz.subvec(evt_gen, smp_gen(evt_idx)), ref_siz, mut_frq_tmp(i), evt_gen, smp_gen(evt_idx), ptn_num);
        mut_frq_pth.submat(i, (evt_gen - smp_gen(0)) * ptn_num, i, (smp_gen(evt_idx) - smp_gen(0)) * ptn_num) = path;
        mut_frq_tmp(i) = arma::as_scalar(path.tail(1));
      }
      gen_frq_tmp = calculateGenoFrq_arma(sel_cof(1), dom_par, mut_frq_tmp);
      wght_tmp = calculateParticleWeight_arma(gen_lik(evt_idx), sel_cof(1), dom_par, mut_frq_tmp);

      if (arma::sum(wght_tmp) > 0) {
        // arma::dcolvec prob = arma::normalise(wght_tmp, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, wght_tmp);

        lik = lik * arma::mean(wght_tmp);
        wght.col(evt_idx) = wght_tmp;
        mut_frq_pre.col(evt_idx) = mut_frq_tmp;
        mut_frq_pst.col(evt_idx) = mut_frq_tmp.elem(indx);
        mut_frq_pth = mut_frq_pth.rows(indx);
        gen_frq_pre.slice(evt_idx) = gen_frq_tmp;
        gen_frq_pst.slice(evt_idx) = gen_frq_tmp.cols(indx);
      } else {
        lik = 0;
        wght.shed_cols(evt_idx, smp_gen.n_elem - 1);
        mut_frq_pre.shed_cols(evt_idx, smp_gen.n_elem - 1);
        mut_frq_pst.shed_cols(evt_idx, smp_gen.n_elem - 1);
        gen_frq_pre.shed_slices(evt_idx, smp_gen.n_elem - 1);
        gen_frq_pst.shed_slices(evt_idx, smp_gen.n_elem - 1);
      }
    } else {
      cout << "generation: " << smp_gen(evt_idx) << endl;
      wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
      mut_frq_tmp = mut_frq_pst.col(evt_idx - 1);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateWFD_arma(sel_cof(1), dom_par, pop_siz.subvec(smp_gen(evt_idx - 1), smp_gen(evt_idx)), ref_siz, mut_frq_tmp(i), smp_gen(evt_idx - 1), smp_gen(evt_idx), ptn_num);
        mut_frq_pth.submat(i, (smp_gen(evt_idx - 1) - smp_gen(0)) * ptn_num, i, (smp_gen(evt_idx) - smp_gen(0)) * ptn_num) = path;
        mut_frq_tmp(i) = arma::as_scalar(path.tail(1));
      }
      gen_frq_tmp = calculateGenoFrq_arma(sel_cof(1), dom_par, mut_frq_tmp);
      wght_tmp = calculateParticleWeight_arma(gen_lik(evt_idx), sel_cof(1), dom_par, mut_frq_tmp);

      if (arma::sum(wght_tmp) > 0) {
        // arma::dcolvec prob = arma::normalise(wght_tmp, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, wght_tmp);

        lik = lik * arma::mean(wght_tmp);
        wght.col(evt_idx) = wght_tmp;
        mut_frq_pre.col(evt_idx) = mut_frq_tmp;
        mut_frq_pst.col(evt_idx) = mut_frq_tmp.elem(indx);
        mut_frq_pth = mut_frq_pth.rows(indx);
        gen_frq_pre.slice(evt_idx) = gen_frq_tmp;
        gen_frq_pst.slice(evt_idx) = gen_frq_tmp.cols(indx);
      } else {
        lik = 0;
        wght.shed_cols(evt_idx, smp_gen.n_elem - 1);
        mut_frq_pre.shed_cols(evt_idx, smp_gen.n_elem - 1);
        mut_frq_pst.shed_cols(evt_idx, smp_gen.n_elem - 1);
        gen_frq_pre.shed_slices(evt_idx, smp_gen.n_elem - 1);
        gen_frq_pst.shed_slices(evt_idx, smp_gen.n_elem - 1);
      }
    }

    for (arma::uword k = evt_idx + 1; k < smp_gen.n_elem; k++) {
      cout << "generation: " << smp_gen(k) << endl;
      wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
      mut_frq_tmp = mut_frq_pst.col(k - 1);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateWFD_arma(sel_cof(1), dom_par, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, mut_frq_tmp(i), smp_gen(k - 1), smp_gen(k), ptn_num);
        mut_frq_pth.submat(i, (smp_gen(k - 1) - smp_gen(0)) * ptn_num, i, (smp_gen(k) - smp_gen(0)) * ptn_num) = path;
        mut_frq_tmp(i) = arma::as_scalar(path.tail(1));
      }
      gen_frq_tmp = calculateGenoFrq_arma(sel_cof(1), dom_par, mut_frq_tmp);
      wght_tmp = calculateParticleWeight_arma(gen_lik(k), sel_cof(1), dom_par, mut_frq_tmp);

      if (arma::sum(wght_tmp) > 0) {
        // arma::dcolvec prob = arma::normalise(wght_tmp, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, wght_tmp);

        lik = lik * arma::mean(wght_tmp);
        wght.col(k) = wght_tmp;
        mut_frq_pre.col(k) = mut_frq_tmp;
        mut_frq_pst.col(k) = mut_frq_tmp.elem(indx);
        mut_frq_pth = mut_frq_pth.rows(indx);
        gen_frq_pre.slice(k) = gen_frq_tmp;
        gen_frq_pst.slice(k) = gen_frq_tmp.cols(indx);
      } else {
        lik = 0;
        wght.shed_cols(k, smp_gen.n_elem - 1);
        mut_frq_pre.shed_cols(k, smp_gen.n_elem - 1);
        mut_frq_pst.shed_cols(k, smp_gen.n_elem - 1);
        gen_frq_pre.shed_slices(k, smp_gen.n_elem - 1);
        gen_frq_pst.shed_slices(k, smp_gen.n_elem - 1);

        break;
      }
    }
  }

  return List::create(Named("lik", lik),
                      Named("wght", wght),
                      Named("mut_frq_pth", mut_frq_pth),
                      Named("mut_frq_pre_resmp", mut_frq_pre),
                      Named("mut_frq_pst_resmp", mut_frq_pst),
                      Named("gen_frq_pre_resmp", gen_frq_pre),
                      Named("gen_frq_pst_resmp", gen_frq_pst));
}
/*************************/


/********** PMMH *********/
// Calculate the log-likelihood using the bootstrap particle filter
// [[Rcpp::export]]
void calculateLogLikelihood_arma(double& log_lik, arma::drowvec& frq_pth, const arma::dcolvec& sel_cof, const double& dom_par, const arma::icolvec& pop_siz, const int& ref_siz, const int& evt_gen, const arma::irowvec& smp_gen, const arma::field<arma::dmat>& gen_lik, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  log_lik = 0;
  frq_pth = arma::zeros<arma::drowvec>(arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1);

  arma::uword evt_idx = arma::as_scalar(arma::max(arma::find(smp_gen <= evt_gen))) + 1;

  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat mut_frq_pth = arma::zeros<arma::dmat>(pcl_num, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1);
  arma::dcolvec mut_frq_pre = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dcolvec mut_frq_pst = arma::zeros<arma::dcolvec>(pcl_num);

  // before the event of interest
  // initialise the particles
  mut_frq_pre = arma::randu<arma::dcolvec>(pcl_num);
  wght = calculateParticleWeight_arma(gen_lik(0), sel_cof(0), dom_par, mut_frq_pre);

  if (arma::mean(wght) > 0) {
    log_lik = log_lik + log(arma::mean(wght));
    // arma::dcolvec prob = arma::normalise(wght, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, wght);
    mut_frq_pst = mut_frq_pre.elem(indx);
  } else {
    log_lik = -(arma::datum::inf);

    return;
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < evt_idx; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::drowvec path = simulateWFD_arma(sel_cof(0), dom_par, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, mut_frq_pst(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      mut_frq_pth.submat(i, (smp_gen(k - 1) - smp_gen(0)) * ptn_num, i, (smp_gen(k) - smp_gen(0)) * ptn_num) = path;
      mut_frq_pre(i) = arma::as_scalar(path.tail(1));
    }
    wght = calculateParticleWeight_arma(gen_lik(k), sel_cof(0), dom_par, mut_frq_pre);

    if (arma::mean(wght) > 0) {
      log_lik = log_lik + log(arma::mean(wght));
      // arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, wght);
      mut_frq_pst = mut_frq_pre.elem(indx);
      mut_frq_pth = mut_frq_pth.rows(indx);
    } else {
      log_lik = -(arma::datum::inf);

      return;
    }
  }

  // after the event of interest
  // run the bootstrap particle filter
  if (evt_idx < smp_gen.n_elem) {
    if (smp_gen(evt_idx - 1) != evt_gen) {
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateWFD_arma(sel_cof(0), dom_par, pop_siz.subvec(smp_gen(evt_idx - 1), evt_gen), ref_siz, mut_frq_pst(i), smp_gen(evt_idx - 1), evt_gen, ptn_num);
        mut_frq_pth.submat(i, (smp_gen(evt_idx - 1) - smp_gen(0)) * ptn_num, i, (evt_gen - smp_gen(0)) * ptn_num) = path;
        mut_frq_pre(i) = arma::as_scalar(path.tail(1));
      }
      mut_frq_pst = mut_frq_pre;

      wght = arma::zeros<arma::dcolvec>(pcl_num);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateWFD_arma(sel_cof(1), dom_par, pop_siz.subvec(evt_gen, smp_gen(evt_idx)), ref_siz, mut_frq_pst(i), evt_gen, smp_gen(evt_idx), ptn_num);
        mut_frq_pth.submat(i, (evt_gen - smp_gen(0)) * ptn_num, i, (smp_gen(evt_idx) - smp_gen(0)) * ptn_num) = path;
        mut_frq_pre(i) = arma::as_scalar(path.tail(1));
      }
      wght = calculateParticleWeight_arma(gen_lik(evt_idx), sel_cof(1), dom_par, mut_frq_pre);

      if (arma::mean(wght) > 0) {
        log_lik = log_lik + log(arma::mean(wght));
        // arma::dcolvec prob = arma::normalise(wght, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, wght);
        mut_frq_pst = mut_frq_pre.elem(indx);
        mut_frq_pth = mut_frq_pth.rows(indx);
      } else {
        log_lik = -(arma::datum::inf);

        return;
      }
    } else {
      wght = arma::zeros<arma::dcolvec>(pcl_num);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateWFD_arma(sel_cof(1), dom_par, pop_siz.subvec(smp_gen(evt_idx - 1), smp_gen(evt_idx)), ref_siz, mut_frq_pst(i), smp_gen(evt_idx - 1), smp_gen(evt_idx), ptn_num);
        mut_frq_pth.submat(i, (smp_gen(evt_idx - 1) - smp_gen(0)) * ptn_num, i, (smp_gen(evt_idx) - smp_gen(0)) * ptn_num) = path;
        mut_frq_pre(i) = arma::as_scalar(path.tail(1));
      }
      wght = calculateParticleWeight_arma(gen_lik(evt_idx), sel_cof(1), dom_par, mut_frq_pre);

      if (arma::mean(wght) > 0) {
        log_lik = log_lik + log(arma::mean(wght));
        // arma::dcolvec prob = arma::normalise(wght, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, wght);
        mut_frq_pst = mut_frq_pre.elem(indx);
        mut_frq_pth = mut_frq_pth.rows(indx);
      } else {
        log_lik = -(arma::datum::inf);

        return;
      }
    }

    for (arma::uword k = evt_idx + 1; k < smp_gen.n_elem; k++) {
      wght = arma::zeros<arma::dcolvec>(pcl_num);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateWFD_arma(sel_cof(1), dom_par, pop_siz.subvec(smp_gen(k - 1), smp_gen(k)), ref_siz, mut_frq_pst(i), smp_gen(k - 1), smp_gen(k), ptn_num);
        mut_frq_pth.submat(i, (smp_gen(k - 1) - smp_gen(0)) * ptn_num, i, (smp_gen(k) - smp_gen(0)) * ptn_num) = path;
        mut_frq_pre(i) = arma::as_scalar(path.tail(1));
      }
      wght = calculateParticleWeight_arma(gen_lik(k), sel_cof(1), dom_par, mut_frq_pre);

      if (arma::mean(wght) > 0) {
        log_lik = log_lik + log(arma::mean(wght));
        // arma::dcolvec prob = arma::normalise(wght, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, wght);
        mut_frq_pst = mut_frq_pre.elem(indx);
        mut_frq_pth = mut_frq_pth.rows(indx);
      } else {
        log_lik = -(arma::datum::inf);

        return;
      }
    }

    frq_pth = mut_frq_pth.row(0);
  } else {
    frq_pth = mut_frq_pth.row(0);
  }

  return;
}

// Calculate the optimal particle number in the particle marginal Metropolis-Hastings
// [[Rcpp::export]]
List calculateOptimalParticleNum_arma(const arma::dcolvec& sel_cof, const double& dom_par, const arma::icolvec& pop_siz, const int& ref_siz, const int& evt_gen, const arma::dmat& raw_smp, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& gap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::irowvec smp_gen = arma::unique(arma::conv_to<arma::irowvec>::from(raw_smp.row(0)));
  arma::field<arma::dmat> gen_lik(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    arma::ucolvec idx = arma::find(raw_smp.row(0) == smp_gen(k));
    gen_lik(k) = raw_smp.submat(1, idx.min(), 3, idx.max());
  }

  arma::drowvec log_lik(300);
  arma::drowvec frq_pth = arma::zeros<arma::drowvec>(arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1);
  for (arma::uword i = 0; i < 300; i++) {
    calculateLogLikelihood_arma(log_lik(i), frq_pth, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, gen_lik, ptn_num, pcl_num);
  }

  arma::drowvec log_lik_sdv(1);
  log_lik_sdv(0) = arma::stddev(log_lik);
  log_lik_sdv.print();
  arma::urowvec opt_pcl_num(1);
  opt_pcl_num(0) = pcl_num;

  if (log_lik_sdv(0) > 1.7) {
    while (log_lik_sdv(0) > 1.0) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) + gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        calculateLogLikelihood_arma(log_lik(i), frq_pth, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, gen_lik, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
    opt_pcl_num = arma::reverse(opt_pcl_num);
    log_lik_sdv = arma::reverse(log_lik_sdv);
  } else if (log_lik_sdv(0) < 1.0) {
    while (log_lik_sdv(0) < 1.7 && opt_pcl_num(0) > gap_num) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) - gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        calculateLogLikelihood_arma(log_lik(i), frq_pth, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, gen_lik, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
  } else {
    while (log_lik_sdv(0) > 1.0) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) + gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        calculateLogLikelihood_arma(log_lik(i), frq_pth, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, gen_lik, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
    opt_pcl_num = arma::reverse(opt_pcl_num);
    log_lik_sdv = arma::reverse(log_lik_sdv);

    while (log_lik_sdv(0) < 1.7 && opt_pcl_num(0) > gap_num) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) - gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        calculateLogLikelihood_arma(log_lik(i), frq_pth, sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, gen_lik, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
  }

  return List::create(Named("opt_pcl_num", opt_pcl_num),
                      Named("log_lik_sdv", log_lik_sdv));
}

// Run the particle marginal Metropolis-Hastings
//[[Rcpp::export]]
List runPMMH_arma(const arma::dcolvec& sel_cof, const double& dom_par, const arma::icolvec& pop_siz, const int& ref_siz, const int& evt_gen, const arma::dmat& raw_smp, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::irowvec smp_gen = arma::unique(arma::conv_to<arma::irowvec>::from(raw_smp.row(0)));
  arma::field<arma::dmat> gen_lik(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    arma::ucolvec idx = arma::find(raw_smp.row(0) == smp_gen(k));
    gen_lik(k) = raw_smp.submat(1, idx.min(), 3, idx.max());
  }

  arma::dmat sel_cof_chn = arma::zeros<arma::dmat>(2, itn_num);
  arma::drowvec dom_par_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::dmat frq_pth_chn = arma::zeros<arma::dmat>(itn_num, arma::uword(arma::max(raw_smp.row(0)) - arma::min(raw_smp.row(0))) + 1);

  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(2);
  arma::drowvec log_lik = arma::zeros<arma::drowvec>(2);
  arma::drowvec frq_pth = arma::zeros<arma::drowvec>(arma::uword(arma::max(raw_smp.row(0)) - arma::min(raw_smp.row(0))) * ptn_num + 1);

  arma::dcolvec sel_cof_sd = {5e-03, 5e-03};
  double dom_par_sd = 5e-03;

  // initialise the samples and population genetic parameters
  cout << "iteration: " << 1 << endl;

  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn.col(0) = sel_cof;
  dom_par_chn(0) = dom_par;

  calculateLogLikelihood_arma(log_lik(0), frq_pth, sel_cof_chn.col(0), dom_par_chn(0), pop_siz, ref_siz, evt_gen, smp_gen, gen_lik, ptn_num, pcl_num);
  arma::ucolvec frq_idx = arma::linspace<arma::ucolvec>(0, smp_gen.max() - smp_gen.min(), smp_gen.max() - smp_gen.min() + 1) * ptn_num;
  frq_pth_chn.row(0) = arma::conv_to<arma::drowvec >::from(frq_pth.elem(frq_idx));

  double apt_cnt = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;

    // draw the candidate of the selection coefficients from the random walk proposal
    sel_cof_chn.col(i) = sel_cof_chn.col(i - 1) + sel_cof_sd % arma::randn<arma::dcolvec>(2);
    dom_par_chn(i) = dom_par_chn(i - 1) + dom_par_sd * arma::randn();

    if (arma::any(sel_cof_chn.col(i) < -1)) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      dom_par_chn(i) = dom_par_chn(i - 1);
      frq_pth_chn.row(i) = frq_pth_chn.row(i - 1);
      // log_lik(1) = log_lik(0);
      // apt_cnt = apt_cnt + 0;
      cout << "acceptance: " << apt_cnt / i << endl;
    // } else if (arma::any(sel_cof_chn.col(i) * dom_par_chn(i) < -1)) {
    } else if (dom_par_chn(i) < 0 || dom_par_chn(i) > 1) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      dom_par_chn(i) = dom_par_chn(i - 1);
      frq_pth_chn.row(i) = frq_pth_chn.row(i - 1);
      // log_lik(1) = log_lik(0);
      // apt_cnt = apt_cnt + 0;
      cout << "acceptance: " << apt_cnt / i << endl;
    } else {
      // calculate the proposal
      // arma::drowvec log_psl = arma::zeros<arma::drowvec>(2);

      // calculate the likelihood
      calculateLogLikelihood_arma(log_lik(1), frq_pth, sel_cof_chn.col(i), dom_par_chn(i), pop_siz, ref_siz, evt_gen, smp_gen, gen_lik, ptn_num, pcl_num);
      frq_pth_chn.row(i) = arma::conv_to<arma::drowvec >::from(frq_pth.elem(frq_idx));

      // calculate the acceptance ratio
      // double apt_rto = exp(log_lik(1) - log_lik(0));
      // double apt_rto = exp((log_pri(1) + log_lik(1) + log_psl(1)) - (log_pri(0) + log_lik(0) + log_psl(0)));

      if (arma::randu() > exp(log_lik(1) - log_lik(0))) {
        sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
        dom_par_chn(i) = dom_par_chn(i - 1);
        frq_pth_chn.row(i) = frq_pth_chn.row(i - 1);
        // log_lik(1) = log_lik(0);
        // apt_cnt = apt_cnt + 0;
        cout << "acceptance: " << apt_cnt / i << endl;
      } else {
        log_lik(0) = log_lik(1);
        apt_cnt = apt_cnt + 1;
        cout << "acceptance: " << apt_cnt / i << endl;
      }
    }
  }

  return List::create(Named("sel_cof_chn", sel_cof_chn),
                      Named("dom_par_chn", dom_par_chn),
                      Named("frq_pth_chn", frq_pth_chn));
}

// Run the adaptive particle marginal Metropolis-Hastings
//[[Rcpp::export]]
List runAdaptPMMH_arma(const arma::dcolvec& sel_cof, const double& dom_par, const arma::icolvec& pop_siz, const int& ref_siz, const int& evt_gen, const arma::dmat& raw_smp, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num, const arma::drowvec& stp_siz, double& apt_rto) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::irowvec smp_gen = arma::unique(arma::conv_to<arma::irowvec>::from(raw_smp.row(0)));
  arma::field<arma::dmat> gen_lik(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    arma::ucolvec idx = arma::find(raw_smp.row(0) == smp_gen(k));
    gen_lik(k) = raw_smp.submat(1, idx.min(), 3, idx.max());
  }

  arma::dmat sel_cof_chn = arma::zeros<arma::dmat>(2, itn_num);
  arma::drowvec dom_par_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::dmat frq_pth_chn = arma::zeros<arma::dmat>(itn_num, arma::uword(arma::max(raw_smp.row(0)) - arma::min(raw_smp.row(0))) + 1);

  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(2);
  arma::drowvec log_lik = arma::zeros<arma::drowvec>(2);
  arma::drowvec frq_pth = arma::zeros<arma::drowvec>(arma::uword(arma::max(raw_smp.row(0)) - arma::min(raw_smp.row(0))) * ptn_num + 1);

  arma::dcolvec U = arma::zeros<arma::dcolvec>(3);
  arma::dcolvec V = arma::zeros<arma::dcolvec>(3);
  arma::dmat S = {{5e-03, 0e-00, 0e-00},
                  {0e-00, 5e-03, 0e-00},
                  {0e-00, 0e-00, 5e-03}};
  arma::dmat M = arma::zeros<arma::dmat>(3, 3);
  arma::dmat I = arma::eye<arma::dmat>(3, 3);

  // initialise the samples and population genetic parameters
  cout << "iteration: " << 1 << endl;

  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn.col(0) = sel_cof;
  dom_par_chn(0) = dom_par;

  calculateLogLikelihood_arma(log_lik(0), frq_pth, sel_cof_chn.col(0), dom_par_chn(0), pop_siz, ref_siz, evt_gen, smp_gen, gen_lik, ptn_num, pcl_num);
  arma::ucolvec frq_idx = arma::linspace<arma::ucolvec>(0, smp_gen.max() - smp_gen.min(), smp_gen.max() - smp_gen.min() + 1) * ptn_num;
  frq_pth_chn.row(0) = arma::conv_to<arma::drowvec >::from(frq_pth.elem(frq_idx));

  double apt_cnt = 0;
  double alpha = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;

    // draw the candidate of the selection coefficients from the random walk proposal
    U = arma::randn<arma::dcolvec>(3);
    V = S * U;
    sel_cof_chn.col(i) = sel_cof_chn.col(i - 1) + V.subvec(0, 1);
    dom_par_chn(i) = dom_par_chn(i - 1) + V(2);

    alpha = 0;
    if (arma::any(sel_cof_chn.col(i) < -1)) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      dom_par_chn(i) = dom_par_chn(i - 1);
      frq_pth_chn.row(i) = frq_pth_chn.row(i - 1);
      // log_lik(1) = log_lik(0);
      // apt_cnt = apt_cnt + 0;
      cout << "acceptance: " << apt_cnt / i << endl;
    // } else if (arma::any(sel_cof_chn.col(i) * dom_par_chn(i) < -1)) {
    } else if (dom_par_chn(i) < 0 || dom_par_chn(i) > 1) {
      sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
      dom_par_chn(i) = dom_par_chn(i - 1);
      frq_pth_chn.row(i) = frq_pth_chn.row(i - 1);
      // log_lik(1) = log_lik(0);
      // apt_cnt = apt_cnt + 0;
      cout << "acceptance: " << apt_cnt / i << endl;
    } else {
      // calculate the proposal
      // arma::drowvec log_psl = arma::zeros<arma::drowvec>(2);

      // calculate the likelihood
      calculateLogLikelihood_arma(log_lik(1), frq_pth, sel_cof_chn.col(i), dom_par_chn(i), pop_siz, ref_siz, evt_gen, smp_gen, gen_lik, ptn_num, pcl_num);
      frq_pth_chn.row(i) = arma::conv_to<arma::drowvec >::from(frq_pth.elem(frq_idx));

      // calculate the acceptance ratio
      // double apt_rto = exp(log_lik(1) - log_lik(0));
      // double apt_rto = exp((log_pri(1) + log_lik(1) + log_psl(1)) - (log_pri(0) + log_lik(0) + log_psl(0)));

      alpha = arma::find_finite(log_lik).is_empty() ? 0 : exp(log_lik(1) - log_lik(0));
      alpha = (alpha > 1) ? 1 : alpha;
      if (arma::randu() > alpha) {
        sel_cof_chn.col(i) = sel_cof_chn.col(i - 1);
        dom_par_chn(i) = dom_par_chn(i - 1);
        frq_pth_chn.row(i) = frq_pth_chn.row(i - 1);
        // log_lik(1) = log_lik(0);
        // apt_cnt = apt_cnt + 0;
        cout << "acceptance: " << apt_cnt / i << endl;
      } else {
        log_lik(0) = log_lik(1);
        apt_cnt = apt_cnt + 1;
        cout << "acceptance: " << apt_cnt / i << endl;
      }
    }

    U = arma::normalise(U);
    M = S * (I + stp_siz(i) * (alpha - apt_rto) * (U * U.t())) * S.t();
    S = arma::chol(M, "lower");
    // cout << M << endl;
  }

  return List::create(Named("sel_cof_chn", sel_cof_chn),
                      Named("dom_par_chn", dom_par_chn),
                      Named("frq_pth_chn", frq_pth_chn));
}
/*************************/
