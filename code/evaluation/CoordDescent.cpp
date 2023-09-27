#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

int sign(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

NumericVector coord_descent(NumericVector prsty, NumericMatrix prs_cov, double alpha, double lambda, NumericVector beta) {
    /* 
    Coordinate descent to find a minimizer.
    Note: after we solve for beta[prs_j], we use its value from then on.
    */
    for (int prs_j = 0; prs_j < prs_cov.nrow(); prs_j++) {
        /* Statistics calculation */
        double A = prsty[prs_j];
        double B = prs_cov(prs_j, prs_j);
        double C = sum((prs_cov(prs_j, _) * beta)) - B * beta[prs_j];
        /*if (alpha != 0 && prs_j == 0) {
            cout << "\nalpha: " << alpha
            << "\nlambda: " << lambda
            << "\nA: " << A
            << "\nB: " << B
            << "\nC: " << C
            << "\nC start: " << sum((prs_cov(prs_j, _) * beta))
            << "\nbeta[prs_j]: " << beta[prs_j];
        }*/
        /* beta update rule */
        if (abs(A - C) <= (alpha * lambda)) {
            beta[prs_j] = 0;
        } else {
            beta[prs_j] = (sign(A - C) * (abs(A - C) - alpha * lambda)) / (B + (1 - alpha) * lambda);
        }
    }

    return beta;
}

// [[Rcpp::export]]
NumericVector wrap_coord_des(NumericVector prsty, NumericVector prs_cov_vector, int prs_cov_size, double alpha, double lambda, NumericVector beta, double thr) {

    /* Convert prs_cov to an Rcpp Matrix */
    prs_cov_vector.attr("dim") = Dimension(prs_cov_size, prs_cov_size);
    NumericMatrix prs_cov = as<NumericMatrix>(prs_cov_vector);

    NumericVector beta_prev = clone(beta);

    /* First iteration to update beta. */
    int ite_num = 0;
    beta = coord_descent(prsty, prs_cov, alpha, lambda, beta);

    /* Update beta iteratively until threshold (thr) is reached. */
    while (max(abs(beta - beta_prev)) > thr) {
        ++ite_num;
        beta_prev = clone(beta);
        beta = coord_descent(prsty, prs_cov, alpha, lambda, beta);
        /* Break out of loop after max 10,000 iterations. */
        if (ite_num == 10000){
            break;
        }
    }

    return(beta);
}