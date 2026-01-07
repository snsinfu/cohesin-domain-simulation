import numpy as np
import scipy.stats


def fit_power_law(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    mask = (x > 0) & (y > 0)
    x = x[mask]
    y = y[mask]
    linreg = scipy.stats.linregress(np.log(x), np.log(y))
    alpha = linreg.slope
    beta = np.exp(linreg.intercept)
    return alpha, beta


class StatisticalLinearRegression:
    """
    Standard scalar-valued OLS with statistical analyses. Attributes:
    - coef: Regression coefficients.
    - coef_stderr: Standard error estimates for the regression coefficients.
    - coef_dof: The degrees of freedom.
    - coef_pvalue: p-values for two-sided t-test for coefficients being zero.
    - mse: Mean-squared error, or the variance estimate of the regression.
    - R2: Coefficient of determination.
    """
    def __init__(self, X, y):
        X = np.asarray(X)
        y = np.asarray(y)
        valid = np.isfinite(X).all(axis=1) & np.isfinite(y)
        X = X[valid]
        y = y[valid]

        XX_inv = np.linalg.pinv(X.T @ X)
        coef = XX_inv @ (X.T @ y)

        mse = np.mean(np.square(y - X @ coef))
        R2 = 1 - mse / y.var()

        coef_cov = mse * XX_inv
        coef_stderr = np.sqrt(np.diag(coef_cov))
        coef_dof = X.shape[0] - X.shape[1]
        coef_tvalue = coef / coef_stderr
        coef_pvalue = (1 - scipy.stats.t.cdf(np.abs(coef_tvalue), coef_dof)) * 2

        self.coef = coef
        self.coef_cov = coef_cov
        self.coef_stderr = coef_stderr
        self.coef_dof = coef_dof
        self.coef_pvalue = coef_pvalue
        self.mse = mse
        self.R2 = R2

    def predict(self, X):
        """
        Predicts response and standard error for given inputs. Returns a
        tuple (responses, standard errors).
        """
        X = np.asarray(X)
        y = X @ self.coef
        y_stderr = np.sqrt(np.diag(X @ self.coef_cov @ X.T))
        return y, y_stderr


def rank_values(xs, fraction=False):
    """
    Transform sequence of n values into their rank values. Integral ranks
    within [0, n) are returned if `fraction` is false, and fractional ranks
    within [0, 1) are returned if `fraction` is true.
    """
    xs = np.asarray(xs)
    rank = np.arange(len(xs))
    rank[np.argsort(xs)] = rank.copy()
    if fraction:
        rank = rank / len(rank)
    return rank
