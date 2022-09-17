\title{
Chapter 13
}

\section{Generalized Method of Moments}

\subsection{Introduction}

One of the most popular estimation methods in applied econometrics is the Generalized Method of Moments (GMM). GMM generalizes classical method of moments by allowing for more equations than unknown parameters (so are overidentified) and by allowing general nonlinear functions of the observations and parameters. Together this allows for a fairly rich and flexible estimation framework. GMM includes as special cases OLS, IV, multivariate regression, and 2SLS. It includes both linear and nonlinear models. In this chapter we focus primarily on linear models.

The GMM label and methods were introduced to econometrics in a seminal paper by Lars Hansen (1982). The ideas and methods build on the work of Amemiya $(1974,1977)$, Gallant (1977), and Gallant and Jorgenson (1979). The ideas are closely related to the contemporeneous work of Halbert White (1980, 1982) and White and Domowitz (1984). The methods are also related to what are called estimating equations in the statistics literature. For a review of the latter see Godambe (1991).

\subsection{Moment Equation Models}

All of the models that have been introduced so far can be written as moment equation models where the population parameters solve a system of moment equations. Moment equation models are broader than the models so far considered and understanding their common structure opens up straightforward techniques to handle new econometric models.

Moment equation models take the following form. Let $g_{i}(\beta)$ be a known $\ell \times 1$ function of the $i^{\text {th }}$ observation and a $k \times 1$ parameter $\beta$. A moment equation model is summarized by the moment equations

$$
\mathbb{E}\left[g_{i}(\beta)\right]=0
$$

and a parameter space $\beta \in B$. For example, in the instrumental variables model $g_{i}(\beta)=Z_{i}\left(Y_{i}-X_{i}^{\prime} \beta\right)$.

In general, we say that a parameter $\beta$ is identified if there is a unique mapping from the data distribution to $\beta$. In the context of the model (13.1) this means that there is a unique $\beta$ satisfying (13.1). Since (13.1) is a system of $\ell$ equations with $k$ unknowns, then it is necessary that $\ell \geq k$ for there to be a unique solution. If $\ell=k$ we say that the model is just identified, meaning that there is just enough information to identify the parameters. If $\ell>k$ we say that the model is overidentified, meaning that there is excess information. If $\ell<k$ we say that the model is underidentified, meaning that there is insufficient information to identify the parameters. In general, we assume that $\ell \geq k$ so the model is either just identified or overidentified. 

\subsection{Method of Moments Estimators}

In this section we consider the just-identified case $\ell=k$.

Define the sample analog of (13.5)

$$
\bar{g}_{n}(\beta)=\frac{1}{n} \sum_{i=1}^{n} g_{i}(\beta) .
$$

The method of moments estimator (MME) $\widehat{\beta}_{\mathrm{mm}}$ is the parameter value which sets $\bar{g}_{n}(\beta)=0$. Thus

$$
\bar{g}_{n}\left(\widehat{\beta}_{\mathrm{mm}}\right)=\frac{1}{n} \sum_{i=1}^{n} g_{i}\left(\widehat{\beta}_{\mathrm{mm}}\right)=0 .
$$

The equations (13.3) are known as the estimating equations as they are the equations which determine the estimator $\widehat{\beta}_{\mathrm{mm}}$.

In some contexts (such as those discussed in the examples below) there is an explicit solution for $\widehat{\beta}_{\mathrm{mm}}$. In other cases the solution must be found numerically.

We now show how most of the estimators discussed so far in the textbook can be written as method of moments estimators.

Mean: Set $g_{i}(\mu)=Y_{i}-\mu$. The MME is $\widehat{\mu}=\frac{1}{n} \sum_{i=1}^{n} Y_{i}$.

Mean and Variance: Set

$$
g_{i}\left(\mu, \sigma^{2}\right)=\left(\begin{array}{c}
Y_{i}-\mu \\
\left(Y_{i}-\mu\right)^{2}-\sigma^{2}
\end{array}\right) .
$$

The MME are $\widehat{\mu}=\frac{1}{n} \sum_{i=1}^{n} Y_{i}$ and $\widehat{\sigma}^{2}=\frac{1}{n} \sum_{i=1}^{n}\left(Y_{i}-\widehat{\mu}\right)^{2}$.

OLS: Set $g_{i}(\beta)=X_{i}\left(Y_{i}-X_{i}^{\prime} \beta\right)$. The MME is $\widehat{\beta}=\left(\boldsymbol{X}^{\prime} \boldsymbol{X}\right)^{-1}\left(\boldsymbol{X}^{\prime} \boldsymbol{Y}\right)$.

OLS and Variance: Set

$$
g_{i}\left(\beta, \sigma^{2}\right)=\left(\begin{array}{c}
X_{i}\left(Y_{i}-X_{i}^{\prime} \beta\right) \\
\left(Y_{i}-X_{i}^{\prime} \beta\right)^{2}-\sigma^{2}
\end{array}\right) \text {. }
$$

The MME is $\widehat{\beta}=\left(\boldsymbol{X}^{\prime} \boldsymbol{X}\right)^{-1}\left(\boldsymbol{X}^{\prime} \boldsymbol{Y}\right)$ and $\widehat{\sigma}^{2}=\frac{1}{n} \sum_{i=1}^{n}\left(Y_{i}-X_{i}^{\prime} \widehat{\beta}\right)^{2}$.

Multivariate Least Squares, vector form: $\operatorname{Set} g_{i}(\beta)=\bar{X}_{i}^{\prime}\left(Y_{i}-\bar{X}_{i} \beta\right)$. The MME is $\widehat{\beta}=\left(\sum_{i=1}^{n} \bar{X}_{i}^{\prime} \bar{X}_{i}\right)^{-1}\left(\sum_{i=1}^{n} \bar{X}_{i} Y_{i}\right)$ which is (11.4).

Multivariate Least Squares, matrix form: Set $g_{i}(\boldsymbol{B})=\operatorname{vec}\left(X_{i}\left(Y_{i}^{\prime}-X_{i}^{\prime} \boldsymbol{B}\right)\right)$. The MME is $\widehat{\boldsymbol{B}}=\left(\sum_{i=1}^{n} X_{i} X_{i}^{\prime}\right)^{-1}\left(\sum_{i=1}^{n} X_{i} Y_{i}^{\prime}\right)$ which is (11.6).

Seemingly Unrelated Regression: Set

$$
g_{i}(\beta, \Sigma)=\left(\begin{array}{c}
\bar{X}_{i} \Sigma^{-1}\left(Y_{i}-\bar{X}_{i}^{\prime} \beta\right) \\
\operatorname{vec}\left(\Sigma-\left(Y_{i}-\bar{X}_{i}^{\prime} \beta\right)\left(Y_{i}-\bar{X}_{i}^{\prime} \beta\right)^{\prime}\right)
\end{array}\right)
$$

The MME is $\widehat{\beta}=\left(\sum_{i=1}^{n} \bar{X}_{i} \widehat{\Sigma}^{-1} \bar{X}_{i}^{\prime}\right)^{-1}\left(\sum_{i=1}^{n} \bar{X}_{i} \widehat{\Sigma}^{-1} Y_{i}\right)$ and $\widehat{\Sigma}=n^{-1} \sum_{i=1}^{n}\left(Y_{i}-\bar{X}_{i}^{\prime} \widehat{\beta}\right)\left(Y_{i}-\bar{X}_{i}^{\prime} \widehat{\beta}\right)^{\prime}$.

IV: Set $g_{i}(\beta)=Z_{i}\left(Y_{i}-X_{i}^{\prime} \beta\right)$. The MME is $\widehat{\beta}=\left(\boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1}\left(\boldsymbol{Z}^{\prime} \boldsymbol{Y}\right)$. Generated Regressors: Set

$$
g_{i}(\beta, \boldsymbol{A})=\left(\begin{array}{c}
\boldsymbol{A}^{\prime} Z_{i}\left(Y_{i}-Z_{i}^{\prime} \boldsymbol{A} \beta\right) \\
\operatorname{vec}\left(Z_{i}\left(X_{i}^{\prime}-Z_{i}^{\prime} \boldsymbol{A}\right)\right)
\end{array}\right)
$$

The MME is $\widehat{\boldsymbol{A}}=\left(\sum_{i=1}^{n} Z_{i} Z_{i}^{\prime}\right)^{-1}\left(\sum_{i=1}^{n} Z_{i} X_{i}^{\prime}\right)$ and $\widehat{\beta}=\left(\widehat{\boldsymbol{A}}^{\prime} \boldsymbol{Z}^{\prime} \boldsymbol{Z} \widehat{\boldsymbol{A}}\right)^{-1}\left(\widehat{\boldsymbol{A}}^{\prime} \boldsymbol{Z}^{\prime} \boldsymbol{Y}\right)$.

A common feature of these examples is that the estimator can be written as the solution to a set of estimating equations (13.3). This provides a common framework which enables a convenient development of a unified distribution theory.

\subsection{Overidentified Moment Equations}

In the instrumental variables model $g_{i}(\beta)=Z_{i}\left(Y_{i}-X_{i}^{\prime} \beta\right)$. Thus (13.2) is

$$
\bar{g}_{n}(\beta)=\frac{1}{n} \sum_{i=1}^{n} g_{i}(\beta)=\frac{1}{n} \sum_{i=1}^{n} Z_{i}\left(Y_{i}-X_{i}^{\prime} \beta\right)=\frac{1}{n}\left(\boldsymbol{Z}^{\prime} \boldsymbol{Y}-\boldsymbol{Z}^{\prime} \boldsymbol{X} \beta\right) .
$$

We have defined the method of moments estimator for $\beta$ as the parameter value which sets $\bar{g}_{n}(\beta)=$ 0 . However, when the model is overidentified (if $\ell>k$ ) this is generally impossible as there are more equations than free parameters. Equivalently, there is no choice of $\beta$ which sets (13.4) to zero. Thus the method of moments estimator is not defined for the overidentified case.

While we cannot find an estimator which sets $\bar{g}_{n}(\beta)$ equal to zero we can try to find an estimator which makes $\bar{g}_{n}(\beta)$ as close to zero as possible.

One way to think about this is to define the vector $\mu=\boldsymbol{Z}^{\prime} \boldsymbol{Y}$, the matrix $\boldsymbol{G}=\boldsymbol{Z}^{\prime} \boldsymbol{X}$ and the "error" $\eta=\mu-\boldsymbol{G} \beta$. Then we can write (13.4) as $\mu=\boldsymbol{G} \beta+\eta$. This looks like a regression equation with the $\ell \times 1$ dependent variable $\mu$, the $\ell \times k$ regressor matrix $\boldsymbol{G}$, and the $\ell \times 1$ error vector $\eta$. The goal is to make the error vector $\eta$ as small as possible. Recalling our knowledge about least squares we deduct that a simple method is to regress $\mu$ on $\boldsymbol{G}$, obtaining $\widehat{\beta}=\left(\boldsymbol{G}^{\prime} \boldsymbol{G}\right)^{-1}\left(\boldsymbol{G}^{\prime} \mu\right)$. This minimizes the sum-of-squares $\eta^{\prime} \eta$. This is certainly one way to make $\eta$ "small".

More generally we know that when errors are non-homogeneous it can be more efficient to estimate by weighted least squares. Thus for some weight matrix $\boldsymbol{W}$ consider the estimator

$$
\widehat{\beta}=\left(\boldsymbol{G}^{\prime} \boldsymbol{W} \boldsymbol{G}\right)^{-1}\left(\boldsymbol{G}^{\prime} \boldsymbol{W} \boldsymbol{\mu}\right)=\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1}\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{Y}\right) .
$$

This minimizes the weighted sum of squares $\eta^{\prime} \boldsymbol{W} \eta$. This solution is known as the generalized method of moments (GMM).

The estimator is typically defined as follows. Given a set of moment equations (13.2) and an $\ell \times \ell$ weight matrix $\boldsymbol{W}>0$ the GMM criterion function is defined as

$$
J(\beta)=n \bar{g}_{n}(\beta)^{\prime} \boldsymbol{W} \bar{g}_{n}(\beta) .
$$

The factor " $n$ " is not important for the definition of the estimator but is convenient for the distribution theory. The criterion $J(\beta)$ is the weighted sum of squared moment equation errors. When $\boldsymbol{W}=\boldsymbol{I}_{\ell}$ then $J(\beta)=n \bar{g}_{n}(\beta)^{\prime} \bar{g}_{n}(\beta)=n\left\|\bar{g}_{n}(\beta)\right\|^{2}$, the square of the Euclidean length. Since we restrict attention to positive definite weight matrices $\boldsymbol{W}$ the criterion $J(\beta)$ is non-negative.

Definition 13.1 The Generalized Method of Moments (GMM) estimator is

$$
\widehat{\beta}_{\mathrm{gmm}}=\underset{\beta}{\operatorname{argmin}} J(\beta) .
$$

Recall that in the just-identified case $k=\ell$ the method of moments estimator $\widehat{\beta}_{\mathrm{mm}}$ solves $\bar{g}_{n}\left(\widehat{\beta}_{\mathrm{mm}}\right)=$ 0 . Hence in this case $J\left(\widehat{\beta}_{\mathrm{mm}}\right)=0$ which means that $\widehat{\beta}_{\mathrm{mm}}$ minimizes $J(\beta)$ and equals $\widehat{\beta}_{\mathrm{gmm}}=\widehat{\beta}_{\mathrm{mm}}$. This means that GMM includes MME as a special case. This implies that all of our results for GMM apply to any method of moments estimator.

In the over-identified case the GMM estimator depends on the choice of weight matrix $\boldsymbol{W}$ and so this is an important focus of the theory. In the just-identified case the GMM estimator simplifies to the MME which does not depend on $\boldsymbol{W}$.

The method and theory of the generalized method of moments was developed in an influential paper by Lars Hansen (1982). This paper introduced the method, its asymptotic distribution, the form of the efficient weight matrix, and tests for overidentification.

\subsection{Linear Moment Models}

One of the great advantages of the moment equation framework is that it allows both linear and nonlinear models. However, when the moment equations are linear in the parameters then we have explicit solutions for the estimates and a straightforward asymptotic distribution theory. Hence we start by confining attention to linear moment equations and return to nonlinear moment equations later. In the examples listed earlier the estimators which have linear moment equations include the sample mean, OLS, multivariate least squares, IV, and 2SLS. The estimates which have nonlinear moment equations include the sample variance, SUR, and generated regressors.

In particular, we focus on the overidentified IV model with moment equations

$$
g_{i}(\beta)=Z_{i}\left(Y_{i}-X_{i}^{\prime} \beta\right)
$$

where $Z_{i}$ is $\ell \times 1$ and $X_{i}$ is $k \times 1$.

\subsection{GMM Estimator}

Given (13.5) the sample moment equations are (13.4). The GMM criterion can be written as

$$
J(\beta)=n\left(\boldsymbol{Z}^{\prime} \boldsymbol{Y}-\boldsymbol{Z}^{\prime} \boldsymbol{X} \beta\right)^{\prime} \boldsymbol{W}\left(\boldsymbol{Z}^{\prime} \boldsymbol{Y}-\boldsymbol{Z}^{\prime} \boldsymbol{X} \beta\right) .
$$

The GMM estimator minimizes $J(\beta)$. The first order conditions are

$$
\begin{aligned}
0 &=\frac{\partial}{\partial \beta} J(\widehat{\beta}) \\
&=2 \frac{\partial}{\partial \beta} \bar{g}_{n}(\widehat{\beta})^{\prime} \boldsymbol{W} \bar{g}_{n}(\widehat{\beta}) \\
&=-2\left(\frac{1}{n} \boldsymbol{X}^{\prime} \boldsymbol{Z}\right) \boldsymbol{W}\left(\frac{1}{n} \boldsymbol{Z}^{\prime}(\boldsymbol{Y}-\boldsymbol{X} \widehat{\beta})\right) .
\end{aligned}
$$

The solution is given as follows.

Theorem 13.1 For the overidentified IV model

$$
\widehat{\beta}_{\mathrm{gmm}}=\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1}\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{Y}\right) .
$$

While the estimator depends on $\boldsymbol{W}$ the dependence is only up to scale. This is because if $\boldsymbol{W}$ is replaced by $c W$ for some $c>0, \widehat{\beta}_{\text {gmm }}$ does not change. When $W$ is fixed by the user we call $\widehat{\beta}_{\text {gmm }}$ ane-step GMM estimator. The formula (13.6) applies for the over-identified $(\ell>k)$ and the just-identified $(\ell=k)$ case. When the model is just-identified then $\boldsymbol{X}^{\prime} \boldsymbol{Z}$ is $k \times k$ so expression (13.6) simplifies to

$$
\widehat{\beta}_{\mathrm{gmm}}=\left(\boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{W}^{-1}\left(\boldsymbol{X}^{\prime} \boldsymbol{Z}\right)^{-1}\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{Y}\right)=\left(\boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1}\left(\boldsymbol{Z}^{\prime} \boldsymbol{Y}\right)=\widehat{\beta}_{\mathrm{iv}}
$$

the IV estimator.

The GMM estimator (13.6) resembles the 2SLS estimator (12.29). In fact they are equal when $\boldsymbol{W}=$ $\left(\boldsymbol{Z}^{\prime} \boldsymbol{Z}\right)^{-1}$. This means that the 2SLS estimator is a one-step GMM estimator for the linear model.

Theorem 13.2 If $\boldsymbol{W}=\left(\boldsymbol{Z}^{\prime} \boldsymbol{Z}\right)^{-1}$ then $\widehat{\beta}_{\mathrm{gmm}}=\widehat{\beta}_{2 \text { sls. }}$ Furthermore, if $k=\ell$ then $\widehat{\beta}_{\mathrm{gmm}}=\widehat{\beta}_{\mathrm{iv}}$

\subsection{Distribution of GMM Estimator}

Let $\boldsymbol{Q}=\mathbb{E}\left[Z X^{\prime}\right]$ and $\Omega=\mathbb{E}\left[Z Z^{\prime} e^{2}\right]$. Then

$$
\left(\frac{1}{n} \boldsymbol{X}^{\prime} \boldsymbol{Z}\right) \boldsymbol{W}\left(\frac{1}{n} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right) \underset{p}{\longrightarrow} \boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}
$$

and

$$
\left(\frac{1}{n} \boldsymbol{X}^{\prime} \boldsymbol{Z}\right) \boldsymbol{W}\left(\frac{1}{\sqrt{n}} \boldsymbol{Z}^{\prime} \boldsymbol{e}\right) \underset{d}{\longrightarrow} \boldsymbol{Q}^{\prime} \boldsymbol{W} \mathrm{N}(0, \Omega)
$$

We conclude:

Theorem 13.3 Asymptotic Distribution of GMM Estimator. Under Assumption 12.2, as $n \rightarrow \infty, \sqrt{n}\left(\widehat{\beta}_{\mathrm{gmm}}-\beta\right) \underset{d}{\longrightarrow} \mathrm{N}\left(0, \boldsymbol{V}_{\beta}\right)$ where

$$
\boldsymbol{V}_{\boldsymbol{\beta}}=\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1}\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \Omega \boldsymbol{W} \boldsymbol{Q}\right)\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1} .
$$

The GMM estimator is asymptotically normal with a "sandwich form" asymptotic variance.

Our derivation treated the weight matrix $W$ as if it is non-random but Theorem $13.3$ applies to the random weight matrix case so long as $\widehat{\boldsymbol{W}}$ converges in probability to a positive definite limit $\boldsymbol{W}$. This may require scaling the weight matrix, for example replacing $\widehat{\boldsymbol{W}}=\left(\boldsymbol{Z}^{\prime} \boldsymbol{Z}\right)^{-1}$ with $\widehat{\boldsymbol{W}}=\left(n^{-1} \boldsymbol{Z}^{\prime} \boldsymbol{Z}\right)^{-1}$. Since rescaling the weight matrix does not affect the estimator this is ignored in implementation. 

\subsection{Efficient GMM}

The asymptotic distribution of the GMM estimator $\widehat{\beta}_{\mathrm{gmm}}$ depends on the weight matrix $\boldsymbol{W}$ through the asymptotic variance $\boldsymbol{V}_{\beta}$. The asymptotically optimal weight matrix $\boldsymbol{W}_{0}$ is that which minimizes $\boldsymbol{V}_{\beta}$. This turns out to be $\boldsymbol{W}_{0}=\Omega^{-1}$. The proof is left to Exercise 13.4.

When the GMM estimator $\widehat{\beta}$ is constructed with $\boldsymbol{W}=\boldsymbol{W}_{0}=\Omega^{-1}$ (or a weight matrix which is a consistent estimator of $\boldsymbol{W}_{0}$ ) we call it the Efficient GMM estimator:

$$
\widehat{\beta}_{\mathrm{gmm}}=\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \Omega^{-1} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1}\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \Omega^{-1} \boldsymbol{Z}^{\prime} \boldsymbol{Y}\right) .
$$

Its asymptotic distribution takes a simpler form than in Theorem 13.3. By substituting $\boldsymbol{W}=\boldsymbol{W}_{0}=\Omega^{-1}$ into (13.7) we find

$$
\boldsymbol{V}_{\beta}=\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \boldsymbol{Q}\right)^{-1}\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \Omega \Omega^{-1} \boldsymbol{Q}\right)\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \boldsymbol{Q}\right)^{-1}=\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \boldsymbol{Q}\right)^{-1} .
$$

This is the asymptotic variance of the efficient GMM estimator.

Theorem 13.4 Asymptotic Distribution of GMM with Efficient Weight Ma-
trix. Under Assumption $12.2$ and $\boldsymbol{W}=\Omega^{-1}$, as $n \rightarrow \infty, \sqrt{n}\left(\widehat{\beta}_{\mathrm{gmm}}-\beta\right) \underset{d}{\mathrm{~d}}$
$\mathrm{~N}\left(0, \boldsymbol{V}_{\beta}\right)$ where $\boldsymbol{V}_{\beta}=\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \boldsymbol{Q}\right)^{-1}$.

Theorem 13.5 Efficient GMM. Under Assumption 12.2, for any $\boldsymbol{W}>0$,

$$
\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1}\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \Omega \boldsymbol{W} \boldsymbol{Q}\right)\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1}-\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \boldsymbol{Q}\right)^{-1} \geq 0
$$

The inequality " $\geq$ " can be replaced with " $>$ " if $W \neq \Omega^{-1}$. Thus if $\widehat{\beta}_{\mathrm{gmm}}$ is the efficient GMM estimator and $\widetilde{\beta}_{\text {gmm }}$ is another GMM estimator, then

$$
\operatorname{avar}\left[\widehat{\beta}_{\mathrm{gmm}}\right] \leq \operatorname{avar}\left[\widetilde{\beta}_{\mathrm{gmm}}\right] .
$$

For a proof see Exercise 13.4.

This means that the smallest possible GMM covariance matrix (in the positive definite sense) is achieved by the efficient GMM weight matrix.

$\boldsymbol{W}_{0}=\Omega^{-1}$ is not known in practice but it can be estimated consistently as we discuss in Section $13.10 .$ For any $\widehat{\boldsymbol{W}} \underset{p}{\rightarrow} \boldsymbol{W}_{0}$ the asymptotic distribution in Theorem $13.4$ is unaffected. Consequently we call any $\widehat{\beta}_{\mathrm{gmm}}$ constructed with an estimate of the efficient weight matrix an efficient GMM estimator.

By "efficient" we mean that this estimator has the smallest asymptotic variance in the class of GMM estimators with this set of moment conditions. This is a weak concept of optimality as we are only considering alternative weight matrices $\widehat{\boldsymbol{W}}$. However, it turns out that the GMM estimator is semiparametrically efficient as shown by Gary Chamberlain (1987). If it is known that $\mathbb{E}\left[g_{i}(\beta)\right]=0$ and this is all that is known this is a semi-parametric problem as the distribution of the data is unknown. Chamberlain showed that in this context no semiparametric estimator (one which is consistent globally for the class of models considered) can have a smaller asymptotic variance than $\left(\boldsymbol{G}^{\prime} \Omega^{-1} \boldsymbol{G}\right)^{-1}$ where $\boldsymbol{G}=\mathbb{E}\left[\frac{\partial}{\partial \beta^{\prime}} g_{i}(\beta)\right]$. Since the GMM estimator has this asymptotic variance it is semiparametrically efficient.

The results in this section show that in the linear model no estimator has better asymptotic efficiency than the efficient linear GMM estimator. No estimator can do better (in this first-order asymptotic sense) without imposing additional assumptions.

\subsection{Efficient GMM versus 2SLS}

For the linear model we introduced 2SLS as a standard estimator for $\beta$. Now we have introduced GMM which includes 2SLS as a special case. Is there a context where 2SLS is efficient?

To answer this question recall that 2SLS is GMM given the weight matrix $\widehat{\boldsymbol{W}}=\left(\boldsymbol{Z}^{\prime} \boldsymbol{Z}\right)^{-1}$ or equivalently $\widehat{\boldsymbol{W}}=\left(n^{-1} \boldsymbol{Z}^{\prime} \boldsymbol{Z}\right)^{-1}$ since scaling doesn't matter. Since $\widehat{\boldsymbol{W}} \underset{p}{\longrightarrow}\left(\mathbb{E}\left[Z Z^{\prime}\right]\right)^{-1}$ this is asymptotically equivalent to the weight matrix $\boldsymbol{W}=\left(\mathbb{E}\left[Z Z^{\prime}\right]\right)^{-1}$. In contrast, the efficient weight matrix takes the form $\left(\mathbb{E}\left[Z Z^{\prime} e^{2}\right]\right)^{-1}$. Now suppose that the structural equation error $e$ is conditionally homoskedastic in the sense that $\mathbb{E}\left[e^{2} \mid Z\right]=\sigma^{2}$. Then the efficient weight matrix equals $\boldsymbol{W}=\left(\mathbb{E}\left[Z Z^{\prime}\right]\right)^{-1} \sigma^{-2}$ or equivalently $W=\left(\mathbb{E}\left[Z Z^{\prime}\right]\right)^{-1}$ since scaling doesn't matter. The latter weight matrix is the same as the 2SLS asymptotic weight matrix. This shows that the 2SLS weight matrix is the efficient weight matrix under conditional homoskedasticity.

Theorem 13.6 Under Assumption $12.2$ and $\mathbb{E}\left[e^{2} \mid Z\right]=\sigma^{2}, \widehat{\beta}_{2 \text { sls }}$ is efficient GMM.

This shows that 2SLS is efficient under homoskedasticity. When homoskedasticity holds there is no reason to use efficient GMM over 2SLS. More broadly, when homoskedasticity is a reasonable approximation then 2SLS will be a reasonable estimator. However, this result also shows that in the general case where the error is conditionally heteroskedastic, 2SLS is inefficient relative to efficient GMM.

\subsection{Estimation of the Efficient Weight Matrix}

To construct the efficient GMM estimator we need a consistent estimator $\widehat{\boldsymbol{W}}$ of $\boldsymbol{W}_{0}=\Omega^{-1}$. The convention is to form an estimator $\widehat{\Omega}$ of $\Omega$ and then set $\widehat{\boldsymbol{W}}=\widehat{\Omega}^{-1}$.

The two-step GMM estimator proceeds by using a one-step consistent estimator of $\beta$ to construct the weight matrix estimator $\widehat{\boldsymbol{W}}$. In the linear model the natural one-step estimator for $\beta$ is 2 SLS. Set $\widetilde{e}_{i}=Y_{i}-X_{i}^{\prime} \widehat{\beta}_{2 s l s}, \widetilde{g}_{i}=g_{i}(\widetilde{\beta})=Z_{i} \widetilde{e}_{i}$, and $\bar{g}_{n}=n^{-1} \sum_{i=1}^{n} \widetilde{g}_{i}$. Two moment estimators of $\Omega$ are

$$
\widehat{\Omega}=\frac{1}{n} \sum_{i=1}^{n} \widetilde{g}_{i} \widetilde{g}_{i}^{\prime}
$$

and

$$
\widehat{\Omega}^{*}=\frac{1}{n} \sum_{i=1}^{n}\left(\widetilde{g}_{i}-\bar{g}_{n}\right)\left(\widetilde{g}_{i}-\bar{g}_{n}\right)^{\prime} .
$$

The estimator (13.8) is an uncentered covariance matrix estimator while the estimator (13.9) is a centered version. Either is consistent when $\mathbb{E}[Z e]=0$ which holds under correct specification. However under misspecification we may have $\mathbb{E}[Z e] \neq 0$. In the latter context $\widehat{\Omega}^{*}$ remains an estimator of var $[Z e]$ while $\widehat{\Omega}$ is an estimator of $\mathbb{E}\left[Z Z^{\prime} e^{2}\right]$. In this sense $\widehat{\Omega}^{*}$ is a robust variance estimator. For some testing problems it turns out to be preferable to use a covariance matrix estimator which is robust to the alternative hypothesis. For these reasons estimator (13.9) is generally preferred. The uncentered estimator (13.8) is more commonly seen in practice since it is the default choice by most packages. It is also worth observing that when the model is just identified then $\bar{g}_{n}=0$ so the two are algebraically identical. The choice of weight matrix may also impact covariance matrix estimation as discussed in Section 13.12.

Given the choice of covariance matrix estimator we set $\widehat{W}=\widehat{\Omega}^{-1}$ or $\widehat{W}=\widehat{\Omega}^{*-1}$. Given this weight matrix we construct the two-step GMM estimator as (13.6) using the weight matrix $\widehat{\boldsymbol{W}}$.

Since the 2SLS estimator is consistent for $\beta$, by arguments nearly identical to those used for covariance matrix estimation we can show that $\widehat{\Omega}$ and $\widehat{\Omega}^{*}$ are consistent for $\Omega$ and thus $\widehat{\boldsymbol{W}}$ is consistent for $\Omega^{-1}$. See Exercise 13.3.

This also means that the two-step GMM estimator satisfies the conditions for Theorem 13.4.

Theorem $13.7$ Under Assumption $12.2$ and $\Omega>0$, if $\widehat{W}=\widehat{\Omega}^{-1}$ or $\widehat{W}=$ $\widehat{\Omega}^{*-1}$ where the latter are defined in (13.8) and (13.9) then as $n \rightarrow \infty$, $\sqrt{n}\left(\widehat{\beta}_{\mathrm{gmm}}-\beta\right) \underset{d}{\longrightarrow} \mathrm{N}\left(0, \boldsymbol{V}_{\beta}\right)$ where $\boldsymbol{V}_{\beta}=\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \boldsymbol{Q}\right)^{-1}$.

This shows that the two-step GMM estimator is asymptotically efficient.

The two-step GMM estimator of the IV regression equation can be computed in Stata using the ivregress gmm command. By default it uses formula (13.8). The centered version (13.9) may be selected using the center option.

\subsection{Iterated GMM}

The asymptotic distribution of the two-step GMM estimator does not depend on the choice of the preliminary one-step estimator. However, the actual value of the estimator depends on this choice and so will the finite sample distribution. This is undesirable and likely inefficient. To remove this dependence we can iterate the estimation sequence. Specifically, given $\widehat{\beta}_{\mathrm{gmm}}$ we can construct an updated weight matrix estimate $\widehat{\boldsymbol{W}}$ and then re-estimate $\widehat{\beta}_{\mathrm{gmm}}$. This updating can be iterated until convergence ${ }^{1}$. The result is called the iterated GMM estimator and is a common implementation of efficient GMM.

Interestingly, B. E. Hansen and Lee (2021) show that the iterated GMM estimator is unaffected if the weight matrix is computed with or without centering. Standard errors and test statistics, however, will be affected by the choice.

The iterated GMM estimator of the IV regression equation can be computed in Stata using the ivregress gmm command using the igmm option.

\subsection{Covariance Matrix Estimation}

An estimator of the asymptotic variance of $\widehat{\beta}_{\mathrm{gmm}}$ can be obtained by replacing the matrices in the asymptotic variance formula by consistent estimators.

${ }^{1}$ In practice, "convergence" obtains when the difference between the estimates at subsequent steps is smaller than a prespecified tolerance. A sufficient condition for convergence is that the sequence is a contraction mapping. Indeed, B. Hansen and Lee (2021) have shown that the iterated GMM estimator generally satisfies this condition in large samples. For the one-step or two-step GMM estimator the covariance matrix estimator is

$$
\widehat{\boldsymbol{V}}_{\beta}=\left(\widehat{\boldsymbol{Q}}^{\prime} \widehat{\boldsymbol{W}} \widehat{\boldsymbol{Q}}\right)^{-1}\left(\widehat{\boldsymbol{Q}}^{\prime} \widehat{\boldsymbol{W}} \widehat{\Omega} \widehat{\boldsymbol{W}} \widehat{\boldsymbol{Q}}\right)\left(\widehat{\boldsymbol{Q}}^{\prime} \widehat{\boldsymbol{W}} \widehat{\boldsymbol{Q}}\right)^{-1}
$$

where $\widehat{\boldsymbol{Q}}=\frac{1}{n} \sum_{i=1}^{n} Z_{i} X_{i}^{\prime}$. The weight matrix is constructed using either the uncentered estimator (13.8) or centered estimator (13.9) with the residuals $\widehat{e}_{i}=Y_{i}-X_{i}^{\prime} \widehat{\beta}_{\mathrm{gmm}}$.

For the efficient iterated GMM estimator the covariance matrix estimator is

$$
\widehat{\boldsymbol{V}}_{\beta}=\left(\widehat{\boldsymbol{Q}}^{\prime} \widehat{\Omega}^{-1} \widehat{\boldsymbol{Q}}\right)^{-1}=\left(\left(\frac{1}{n} \boldsymbol{X}^{\prime} \boldsymbol{Z}\right) \widehat{\Omega}^{-1}\left(\frac{1}{n} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)\right)^{-1} .
$$

$\widehat{\Omega}$ can be computed using either the uncentered estimator (13.8) or centered estimator (13.9). Based on the asymptotic approximation the estimator (13.11) can be used as well for the two-step estimator but should use the final residuals $\widehat{e}_{i}=Y_{i}-X_{i}^{\prime} \widehat{\beta}_{\mathrm{gmm}}$.

Asymptotic standard errors are given by the square roots of the diagonal elements of $n^{-1} \widehat{\boldsymbol{V}}_{\beta}$.

It is unclear if it is preferred to use the covariance matrix estimator based on the centered or uncentered estimator of $\Omega$ to construct the covariance matrix estimator. Using the centered estimator results in a smaller covariance matrix and standard errors and thus more "significant" tests based on asymptotic critical values. In contrast the uncentered estimator of $\Omega$ will result in larger standard errors and will thus be more "conservative".

In Stata, the default covariance matrix estimation method is determined by the choice of weight matrix. Thus if the centered estimator (13.9) is used for the weight matrix it is also used for the covariance matrix estimator.

\subsection{Clustered Dependence}

In Section $4.21$ we introduced clustered dependence and in Section $12.25$ described covariance matrix estimation for 2SLS. The methods extend naturally to GMM but with the additional complication of potentially altering weight matrix calculation.

The structural equation for the $g^{t h}$ cluster can be written as the matrix system $\boldsymbol{Y}_{g}=\boldsymbol{X}_{g} \beta+\boldsymbol{e}_{g}$. Using this notation the centered GMM estimator with weight matrix $\boldsymbol{W}$ can be written as

$$
\widehat{\beta}_{\mathrm{gmm}}-\beta=\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W}\left(\sum_{g=1}^{G} \boldsymbol{Z}_{g}^{\prime} \boldsymbol{e}_{g}\right)
$$

The cluster-robust covariance matrix estimator for $\widehat{\beta}_{\mathrm{gmm}}$ is

$$
\widehat{\boldsymbol{V}}_{\beta}=\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \widehat{\boldsymbol{S}} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1}
$$

with

$$
\widehat{\boldsymbol{S}}=\sum_{g=1}^{G} \boldsymbol{Z}_{g}^{\prime} \widehat{\boldsymbol{e}}_{g} \widehat{\boldsymbol{e}}_{g}^{\prime} \boldsymbol{Z}_{g}
$$

and the clustered residuals

$$
\widehat{\boldsymbol{e}}_{g}=\boldsymbol{Y}_{g}-\boldsymbol{X}_{g} \widehat{\beta}_{\mathrm{gmm}} .
$$

The cluster-robust estimator (13.12) is appropriate for the one-step or two-step GMM estimator. It is also appropriate for the iterated estimator when the latter uses a conventional (non-clustered) efficient weight matrix. However in the clustering context it is more natural to use a cluster-robust weight matrix such as $\boldsymbol{W}=\widehat{\boldsymbol{S}}^{-1}$ where $\widehat{\boldsymbol{S}}$ is a cluster-robust covariance estimator as in (13.13) based on a one-step or iterated residual. This gives rise to the cluster-robust GMM estimator

$$
\widehat{\beta}_{\mathrm{gmm}}=\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \widehat{\boldsymbol{S}}^{-1} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{X}^{\prime} \boldsymbol{Z} \widehat{\boldsymbol{S}}^{-1} \boldsymbol{Z}^{\prime} \boldsymbol{Y}
$$

An appropriate cluster-robust covariance matrix estimator is

$$
\widehat{\boldsymbol{V}}_{\beta}=\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \widehat{\boldsymbol{S}}^{-1} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1}
$$

where $\widehat{S}$ is calculated using the final residuals.

To implement a cluster-robust weight matrix use the 2SLS estimator for first step. Compute the cluster residuals (13.14) and covariance matrix (13.13). Then (13.15) is the two-step GMM estimator. Iterating the residuals and covariance matrix until convergence we obtain the iterated GMM estimator.

In Stata, using the ivregress gmm command with the cluster option implements the two-step GMM estimator using the cluster-robust weight matrix and cluster-robust covariance matrix estimator. To use the centered covariance matrix use the center option and to implement the iterated GMM estimator use the igmm option. Alternatively, you can use the wmatrix and vce options to separately specify the weight matrix and covariance matrix estimation methods.

\subsection{Wald Test}

For a given function $r(\beta): \mathbb{R}^{k} \rightarrow \Theta \subset \mathbb{R}^{q}$ we define the parameter $\theta=r(\beta)$. The GMM estimator of $\theta$ is $\widehat{\theta}_{\mathrm{gmm}}=r\left(\widehat{\beta}_{\mathrm{gmm}}\right)$. By the delta method it is asymptotically normal with covariance matrix $\boldsymbol{V}_{\theta}=$ $\boldsymbol{R}^{\prime} \boldsymbol{V}_{\beta} \boldsymbol{R}$ where $\boldsymbol{R}=\frac{\partial}{\partial \beta} r(\beta)^{\prime}$. An estimator of the asymptotic covariance matrix is $\widehat{\boldsymbol{V}}_{\theta}=\widehat{\boldsymbol{R}}^{\prime} \widehat{\boldsymbol{V}}_{\beta} \widehat{\boldsymbol{R}}^{\text {where }}$ $\widehat{\boldsymbol{R}}=\frac{\partial}{\partial \beta} r\left(\widehat{\beta}_{\mathrm{gmm}}\right)^{\prime}$. When $\theta$ is scalar then an asymptotic standard error for $\widehat{\theta}_{\mathrm{gmm}}$ is formed as $\sqrt{n^{-1} \widehat{\boldsymbol{V}}_{\theta}} .$

A standard test of the hypothesis $\mathbb{M}_{0}: \theta=\theta_{0}$ against $\mathbb{M}_{1}: \theta \neq \theta_{0}$ is based on the Wald statistic

$$
W=n\left(\widehat{\theta}-\theta_{0}\right)^{\prime} \widehat{\boldsymbol{V}}_{\widehat{\theta}}^{-1}\left(\widehat{\theta}-\theta_{0}\right) .
$$

Let $G_{q}(u)$ denote the $\chi_{q}^{2}$ distribution function.

Theorem $13.8$ Under Assumption 12.2, Assumption 7.3, and $\mathbb{H}_{0}$, as $n \rightarrow \infty$, $W \underset{d}{\longrightarrow} \chi_{q}^{2}$. For $c$ satisfying $\alpha=1-G_{q}(c), \mathbb{P}\left[W>c \mid \mathbb{H}_{0}\right] \longrightarrow \alpha$ so the test "Reject $\mathbb{H}_{0}$ if $W>c$ " has asymptotic size $\alpha$.

For a proof see Exercise 13.5.

In Stata, the commands test and testparm can be used after ivregress gmm to implement Wald tests of linear hypotheses. The commands nlcom and testnl can be used after ivregress gmm to implement Wald tests of nonlinear hypotheses. 

\subsection{Restricted GMM}

It is often desirable to impose restrictions on the coefficients. In this section we consider estimation subject to the linear constraints $\boldsymbol{R}^{\prime} \beta=\boldsymbol{c}$. In the following section we consider nonlinear constraints.

The constrained GMM estimator minimizes the GMM criterion subject to the constraint. It is

$$
\widehat{\beta}_{\mathrm{cgmm}}=\underset{\boldsymbol{R}^{\prime} \beta=\boldsymbol{c}}{\operatorname{argmin}} J(\beta) .
$$

This is the parameter vector which makes the estimating equations as close to zero as possible with respect to the weighted quadratic distance while imposing the restriction on the parameters.

Suppose the weight matrix $\boldsymbol{W}$ is fixed. Using the methods of Chapter 8 it is straightforward to derive that the constrained GMM estimator is

$$
\widehat{\beta}_{\mathrm{cgmm}}=\widehat{\beta}_{\mathrm{gmm}}-\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{R}\left(\boldsymbol{R}^{\prime}\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{R}\right)^{-1}\left(\boldsymbol{R}^{\prime} \widehat{\beta}_{\mathrm{gmm}}-\boldsymbol{c}\right) .
$$

(For details, see Exercise 13.6.)

We derive the asymptotic distribution under the assumption that the restriction is true. Make the substitution $\boldsymbol{c}=\boldsymbol{R}^{\prime} \beta$ in (13.16) and reorganize to find

$$
\sqrt{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\beta\right)=\left(\boldsymbol{I}_{k}-\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{R}\left(\boldsymbol{R}^{\prime}\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime}\right) \sqrt{n}\left(\widehat{\beta}_{\mathrm{gmm}}-\beta\right)
$$

This is a linear function of $\sqrt{n}\left(\widehat{\beta}_{\mathrm{gmm}}-\beta\right)$. Since the asymptotic distribution of the latter is known the asymptotic distribution of $\sqrt{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\beta\right)$ is a linear function of the former.

$$
\begin{aligned}
&\text { Theorem 13.9 Under Assumptions } 12.2 \text { and 8.3, for the constrained GMM es- } \\
&\text { timator (13.16), } \sqrt{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\beta\right) \underset{d}{\longrightarrow} \mathrm{N}\left(0, \boldsymbol{V}_{\mathrm{cgmm}}\right) \text { as } n \rightarrow \infty \text {, where } \\
&\boldsymbol{V}_{\mathrm{cgmm}}=\boldsymbol{V}_{\beta}-\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1} \boldsymbol{R}\left(\boldsymbol{R}^{\prime}\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime} \boldsymbol{V}_{\beta} \\
&\quad-\boldsymbol{V}_{\beta} \boldsymbol{R}\left(\boldsymbol{R}^{\prime}\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime}\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1} \\
& \\
&+\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1} \boldsymbol{R}\left(\boldsymbol{R}^{\prime}\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime} \boldsymbol{V}_{\beta} \boldsymbol{R}\left(\boldsymbol{R}^{\prime}\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime}\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1}
\end{aligned}
$$

For a proof, see Exercise 13.8. Unfortunately the asymptotic covariance matrix formula (13.18) is quite tedious!

Now suppose that the the weight matrix is set as $W=\widehat{\Omega}^{-1}$, the efficient weight matrix from unconstrained estimation. In this case the constrained GMM estimator can be written as

$$
\widehat{\beta}_{\mathrm{cgmm}}=\widehat{\beta}_{\mathrm{gmm}}-\widehat{\boldsymbol{V}}_{\beta} \boldsymbol{R}\left(\boldsymbol{R}^{\prime} \widehat{\boldsymbol{V}}_{\boldsymbol{\beta}} \boldsymbol{R}\right)^{-1}\left(\boldsymbol{R}^{\prime} \widehat{\beta}_{\mathrm{gmm}}-\boldsymbol{c}\right)
$$

which is the same formula (8.25) as efficient minimum distance. (For details, see Exercise 13.7.) We find that the asymptotic covariance matrix simplifies considerably. Theorem 13.10 Under Assumptions $12.2$ and 8.3, for the efficient constrained GMM estimator (13.19), $\sqrt{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\beta\right) \underset{d}{\longrightarrow} \mathrm{N}\left(0, \boldsymbol{V}_{\mathrm{cgmm}}\right)$ as $n \rightarrow \infty$, where

$$
\boldsymbol{V}_{\mathrm{cgmm}}=\boldsymbol{V}_{\beta}-\boldsymbol{V}_{\beta} \boldsymbol{R}\left(\boldsymbol{R}^{\prime} \boldsymbol{V}_{\beta} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime} \boldsymbol{V}_{\beta} .
$$

For a proof, see Exercise 13.9.

The asymptotic covariance matrix (13.20) can be estimated by

$$
\begin{aligned}
\widehat{\boldsymbol{V}}_{\mathrm{cgmm}} &=\widetilde{\boldsymbol{V}}_{\beta}-\widetilde{\boldsymbol{V}}_{\beta} \boldsymbol{R}\left(\boldsymbol{R}^{\prime} \widetilde{\boldsymbol{V}}_{\beta} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime} \widetilde{\boldsymbol{V}}_{\beta} . \\
\widetilde{\boldsymbol{V}}_{\beta} &=\left(\widehat{\boldsymbol{Q}}^{\prime} \widetilde{\Omega}^{-1} \widehat{\boldsymbol{Q}}\right)^{-1} \\
\widetilde{\Omega} &=\frac{1}{n} \sum_{i=1}^{n} Z_{i} Z_{i}^{\prime} \widetilde{e}_{i}^{2} \\
\widetilde{e}_{i} &=Y_{i}-X_{i}^{\prime} \widehat{\beta}_{\mathrm{cgmm}} .
\end{aligned}
$$

The covariance matrix (13.18) can be estimated similarly, though using (13.10) to estimate $\boldsymbol{V}_{\beta}$. The covariance matrix estimator $\widetilde{\Omega}$ can also be replaced with a centered version.

A constrained iterated GMM estimator can be implemented by setting $\boldsymbol{W}=\widetilde{\Omega}^{-1}$ where $\widetilde{\Omega}$ is defined in (13.22) and then iterating until convergence. This is a natural estimator as it is the appropriate implementation of iterated GMM.

Since both $\widehat{\Omega}$ and $\widetilde{\Omega}$ converge to the same limit $\Omega$ under the assumption that the constraint is true the constrained iterated GMM estimator has the asymptotic distribution given in Theorem 13.10.

\subsection{Nonlinear Restricted GMM}

Nonlinear constraints on the parameters can be written as $r(\beta)=0$ for some function $r: \mathbb{R}^{k} \rightarrow \mathbb{R}^{q}$. The constraint is nonlinear if $r(\beta)$ cannot be written as a linear function of $\beta$. Least squares estimation subject to nonlinear constraints was explored in Section 8.14. In this section we introduce GMM estimation subject to nonlinear constraints.

The constrained GMM estimator minimizes the GMM criterion subject to the constraint. It is

$$
\widehat{\beta}_{\mathrm{cgmm}}=\underset{r(\beta)=0}{\operatorname{argmin}} J(\beta) .
$$

This is the parameter vector which makes the estimating equations as close to zero as possible with respect to the weighted quadratic distance while imposing the restriction on the parameters.

In general there is no explicit solution for $\widehat{\beta}_{\mathrm{cgmm}}$. Instead the solution is found numerically. Fortunately there are excellent nonlinear constrained optimization solvers implemented in standard software packages.

For the asymptotic distribution assume that the restriction $r(\beta)=0$ is true. Using the same methods as in the proof of Theorem $8.10$ we can show that (13.17) approximately holds in the sense that

$$
\sqrt{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\beta\right)=\left(\boldsymbol{I}_{k}-\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{R}\left(\boldsymbol{R}^{\prime}\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{W} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime}\right) \sqrt{n}\left(\widehat{\beta}_{\mathrm{gmm}}-\beta\right)+o_{p}(1)
$$

where $\boldsymbol{R}=\frac{\partial}{\partial \beta} r(\beta)^{\prime}$. Thus the asymptotic distribution of the constrained estimator takes the same form as in the linear case. Theorem $13.11$ Under Assumptions $12.2$ and 8.3, for the constrained GMM estimator (13.23), $\sqrt{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\beta\right) \underset{d}{\longrightarrow} \mathrm{N}\left(0, V_{\mathrm{cgmm}}\right)$ as $n \rightarrow \infty$, where $\boldsymbol{V}_{\mathrm{cgmm}}$ equals (13.18). If $W=\widehat{\Omega}^{-1}$, then $V_{\text {cgmm }}$ equals (13.20).

The asymptotic covariance matrix in the efficient case is estimated by (13.21) with $\boldsymbol{R}$ replaced with $\widehat{\boldsymbol{R}}=\frac{\partial}{\partial \beta} r\left(\widehat{\beta}_{\text {cgmm }}\right)^{\prime}$. The asymptotic covariance matrix (13.18) in the general case is estimated similarly.

To implement an iterated restricted GMM estimator the weight matrix may be set as $\boldsymbol{W}=\widetilde{\Omega}^{-1}$ where $\widetilde{\Omega}$ is defined in (13.22), and then iterated until convergence.

\subsection{Constrained Regression}

Take the conventional projection model $Y=X^{\prime} \beta+e$ with $\mathbb{E}[X e]=0$. This is a special case of GMM as it is model (13.5) with $Z=X$. The just-identified GMM estimator equals least squares $\widehat{\beta}_{\mathrm{gmm}}=\widehat{\beta}_{\mathrm{ols}}$.

In Chapter 8 we discussed estimation of the projection model subject to linear constraints $\boldsymbol{R}^{\prime} \beta=\boldsymbol{c}$, which includes exclusion restrictions. Since the projection model is a special case of GMM the constrained projection model is also constrained GMM. From the results of Section $13.15$ we find that the efficient constrained GMM estimator is

$$
\widehat{\beta}_{\mathrm{cgmm}}=\widehat{\beta}_{\mathrm{ols}}-\widehat{\boldsymbol{V}}_{\beta} \boldsymbol{R}\left(\boldsymbol{R}^{\prime} \widehat{\boldsymbol{V}}_{\beta} \boldsymbol{R}\right)^{-1}\left(\boldsymbol{R}^{\prime} \widehat{\beta}_{\mathrm{ols}}-\boldsymbol{c}\right)=\widehat{\beta}_{\mathrm{emd}},
$$

the efficient minimum distance estimator. Thus for linear constraints on the linear projection model efficient GMM equals efficient minimum distance. Thus one convenient method to implement efficient minimum distance is GMM.

\subsection{Multivariate Regression}

GMM methods can simplify estimation and inference for multivariate regressions such as those introduced in Chapter $11 .$

The general multivariate regression (projection) model is

$$
\begin{aligned}
Y_{j} &=X_{j}^{\prime} \beta_{j}+e_{j} \\
\mathbb{E}\left[X_{j} e_{j}\right] &=0
\end{aligned}
$$

for $j=1, \ldots, m$. Using the notation from Section $11.2$ the equations can be written jointly as $Y=\bar{X} \beta+e$ and for the full sample as $\boldsymbol{Y}=\overline{\boldsymbol{X}} \beta+\boldsymbol{e}$. The $\bar{k}$ moment conditions are

$$
\mathbb{E}\left[\bar{X}^{\prime}(Y-\bar{X} \beta)\right]=0 .
$$

Given a $\bar{k} \times \bar{k}$ weight matrix $\boldsymbol{W}$ the GMM criterion is

$$
J(\beta)=n(\boldsymbol{Y}-\overline{\boldsymbol{X}} \beta)^{\prime} \overline{\boldsymbol{X}} \boldsymbol{W} \overline{\boldsymbol{X}}^{\prime}(\boldsymbol{Y}-\overline{\boldsymbol{X}} \beta) .
$$

The GMM estimator $\widehat{\beta}_{\text {gmm }}$ minimizes $J(\beta)$. Since this is a just-identified model the estimator solves the sample equations

$$
\overline{\boldsymbol{X}}^{\prime}\left(\boldsymbol{Y}-\overline{\boldsymbol{X}} \widehat{\beta}_{\mathrm{gmm}}\right)=0
$$

The solution is

$$
\widehat{\beta}_{\mathrm{gmm}}=\left(\sum_{i=1}^{n} \bar{X}_{i}^{\prime} \bar{X}_{i}\right)^{-1}\left(\sum_{i=1}^{n} \bar{X}_{i}^{\prime} Y_{i}\right)=\left(\overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}}\right)^{-1}\left(\overline{\boldsymbol{X}}^{\prime} \boldsymbol{Y}\right)=\widehat{\beta}_{\mathrm{ols}}
$$

the multivariate least squares estimator.

Thus the unconstrained GMM estimator of the multivariate regression model is least squares. The estimator does not depend on the weight matrix since the model is just-identified.

A important advantage of the GMM framework is the ability to incorporate cross-equation constraints. Consider the class of restrictions $\boldsymbol{R}^{\prime} \beta=\boldsymbol{c}$. Minimization of the GMM criterion subject to this restriction has solutions as described in (13.15). The restricted GMM estimator is

$$
\widehat{\beta}_{\mathrm{gmm}}=\widehat{\beta}_{\mathrm{ols}}-\left(\overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}} \boldsymbol{W} \overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}}\right)^{-1} \boldsymbol{R}\left(\boldsymbol{R}^{\prime}\left(\overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}} \boldsymbol{W} \overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}}\right)^{-1} \boldsymbol{R}\right)^{-1}\left(\boldsymbol{R}^{\prime} \widehat{\beta}_{\mathrm{ols}}-\boldsymbol{c}\right)
$$

This estimator depends on the weight matrix because it is over-identified.

A simple choice for weight matrix is $\boldsymbol{W}=\overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}}$. This leads to the one-step estimator

$$
\widehat{\beta}_{1}=\widehat{\beta}_{\mathrm{ols}}-\left(\overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}}\right)^{-1} \boldsymbol{R}\left(\boldsymbol{R}^{\prime}\left(\overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}}\right)^{-1} \boldsymbol{R}\right)^{-1}\left(\boldsymbol{R}^{\prime} \widehat{\beta}_{\mathrm{ols}}-\boldsymbol{c}\right) .
$$

The asymptotically efficient choice sets $\boldsymbol{W}=\widehat{\Omega}^{-1}$ where $\widehat{\Omega}=n^{-1} \sum_{i=1}^{n} \bar{X}_{i}^{\prime} \widehat{e}_{i} \widehat{e}_{i}^{\prime} \bar{X}_{i}$ and $\widehat{e}_{i}=Y_{i}-\bar{X}_{i} \widehat{\beta}_{1}$. This leads to the two-step estimator

$$
\widehat{\beta}_{2}=\widehat{\beta}_{\text {ols }}-\left(\overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}} \widehat{\Omega}^{-1} \overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}}\right)^{-1} \boldsymbol{R}\left(\boldsymbol{R}^{\prime}\left(\overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}} \widehat{\Omega}^{-1} \overline{\boldsymbol{X}}^{\prime} \overline{\boldsymbol{X}}\right)^{-1} \boldsymbol{R}\right)^{-1}\left(\boldsymbol{R}^{\prime} \widehat{\beta}_{\mathrm{ols}}-\boldsymbol{c}\right)
$$

When the regressors $X$ are common across all equations the multivariate regression model can be written conveniently as in (11.3): $Y=\boldsymbol{B}^{\prime} X+e$ with $\mathbb{E}\left[X e^{\prime}\right]=0$. The moment restrictions can be written as the matrix system $\mathbb{E}\left[X\left(Y^{\prime}-X^{\prime} \boldsymbol{B}\right)\right]=0$. Written as a vector system this is (13.24) and leads to the same restricted GMM estimators.

These are general formula for imposing restrictions. In specific cases (such as an exclusion restriction) direct methods may be more convenient. In all cases the solution is found by minimization of the GMM criterion $J(\beta)$ subject to the restriction.

\subsection{Distance Test}

In Section $13.14$ we introduced Wald tests of the hypothesis $\mathbb{M}_{0}: \theta=\theta_{0}$ where $\theta=r(\beta)$ for a given function $r(\beta): \mathbb{R}^{k} \rightarrow \Theta \subset \mathbb{R}^{q}$. When $r(\beta)$ is nonlinear an alternative is to use a criterion-based statistic. This is sometimes called the GMM Distance statistic and sometimes called a LR-like statistic (the LR is for likelihood-ratio). The idea was first put forward by Newey and West (1987a).

The idea is to compare the unrestricted and restricted estimators by contrasting the criterion functions. The unrestricted estimator takes the form

$$
\widehat{\beta}_{\mathrm{gmm}}=\underset{\beta}{\operatorname{argmin}} \widehat{J}(\beta)
$$

where

$$
\widehat{J}(\beta)=n \bar{g}_{n}(\beta)^{\prime} \widehat{\Omega}^{-1} \bar{g}_{n}(\beta)
$$

is the unrestricted GMM criterion with an efficient weight matrix estimate $\widehat{\Omega}$. The minimized value of the criterion is $\widehat{J}=\widehat{J}\left(\widehat{\beta}_{\mathrm{gmm}}\right)$. As in Section 13.15, the estimator subject to $r(\beta)=\theta_{0}$ is

$$
\widehat{\beta}_{\mathrm{cgmm}}=\underset{r(\beta)=\theta_{0}}{\operatorname{argmin}} \widetilde{J}(\beta)
$$

where

$$
\widetilde{J}(\beta)=n \bar{g}_{n}(\beta)^{\prime} \widetilde{\Omega}^{-1} \bar{g}_{n}(\beta)
$$

which depends on an efficient weight matrix estimator, either $\widehat{\Omega}$ (the same as the unrestricted estimator) or $\widetilde{\Omega}$ (the iterated weight matrix from constrained estimation). The minimized value of the criterion is $\widetilde{J}=\widetilde{J}\left(\widehat{\beta}_{\operatorname{cgmm}}\right)$

The GMM distance (or LR-like) statistic is the difference in the criterion functions: $D=\widetilde{J}-\widehat{J}$. The distance test shares the useful feature of LR tests in that it is a natural by-product of the computation of alternative models.

The test has the following large sample distribution.

Theorem $13.12$ Under Assumption 12.2, Assumption 7.3, and $\mathbb{H}_{0}$, then as $n \rightarrow$ $\infty, D \longrightarrow \chi_{q}^{2}$. For $c$ satisfying $\alpha=1-G_{q}(c), \mathbb{P}\left[D>c \mid \mathbb{H}_{0}\right] \longrightarrow \alpha$. The test "Reject $\mathbb{H}_{0}$ if $D>c$ " has asymptotic size $\alpha$.

The proof is given in Section 13.28.

Theorem $13.12$ shows that the distance statistic has the same asymptotic distribution as Wald and likelihood ratio statistics and can be interpreted similarly. Small values of $D$ mean that imposing the restriction does not result in a large value of the moment equations. Hence the restriction appears to be compatible with the data. On the other hand, large values of $D$ mean that imposing the restriction results in a much larger value of the moment equations, implying that the restriction is not compatible with the data. The finding that the asymptotic distribution is chi-squared allows the calculation of asymptotic critical values and p-values.

We now discuss the choice of weight matrix. As mentioned above one simple choice is to set $\widetilde{\Omega}=\widehat{\Omega}$. In this case we have the following result.

Theorem 13.13 If $\widetilde{\Omega}=\widehat{\Omega}$ then $D \geq 0$. Furthermore, if $r$ is linear in $\beta$ then $D$ equals the Wald statistic.

The statement that $\widetilde{\Omega}=\widehat{\Omega}$ implies $D \geq 0$ follows from the fact that in this case the criterion functions $\widehat{J}(\beta)=\widetilde{J}(\beta)$ are identical so the constrained minimum cannot be smaller than the unconstrained. The statement that linear hypotheses and $\widetilde{\Omega}=\widehat{\Omega}$ implies $D=W$ follows from applying the expression for the constrained GMM estimator (13.19) and using the covariance matrix formula (13.11).

The fact that $D \geq 0$ when $\widetilde{\Omega}=\widehat{\Omega}$ motivated Newey and West (1987a) to recommend this choice. However, $\widetilde{\Omega}=\widehat{\Omega}$ is not necessary. Instead, setting $\widetilde{\Omega}$ to equal the constrained efficient weight matrix is natural for efficient estimation of $\widehat{\beta}_{\mathrm{cgmm}}$. In the event that $D<0$ the test simply fails to reject $\mathbb{H}_{0}$ at any significance level.

As discussed in Section $9.17$ for tests of nonlinear hypotheses the Wald statistic can work quite poorly. In particular, the Wald statistic is affected by how the hypothesis $r(\beta)$ is formulated. In contrast, the distance statistic $D$ is not affected by the algebraic formulation of the hypothesis. Current evidence suggests that the $D$ statistic appears to have good sampling properties, and is a preferred test statistic relative to the Wald statistic for nonlinear hypotheses. (See B. E. Hansen (2006).)

In Stata the command estat overid after ivregress gmm can be used to report the value of the GMM criterion $J$. By estimating the two nested GMM regressions the values $\widehat{J}$ and $\widetilde{J}$ can be obtained and $D$ computed.

\subsection{Continuously-Updated GMM}

An alternative to the two-step GMM estimator can be constructed by letting the weight matrix be an explicit function of $\beta$. These leads to the criterion function

$$
J(\beta)=n \bar{g}_{n}(\beta)^{\prime}\left(\frac{1}{n} \sum_{i=1}^{n} g_{i}(\beta) g_{i}(\beta)^{\prime}\right)^{-1} \bar{g}_{n}(\beta) .
$$

The $\widehat{\beta}$ which minimizes this function is called the continuously-updated GMM (CU-GMM) estimator and was introduced by L. Hansen, Heaton and Yaron (1996).

A complication is that the continuously-updated criterion $J(\beta)$ is not quadratic in $\beta$. This means that minimization requires numerical methods. It may appear that the CU-GMM estimator is the same as the iterated GMM estimator but this is not the case at all. They solve distinct first-order conditions and can be quite different in applications.

Relative to traditional GMM the CU-GMM estimator has lower bias but thicker distributional tails. While it has received considerable theoretical attention it is not used commonly in applications.

\subsection{OverIdentification Test}

In Section $12.31$ we introduced the Sargan (1958) overidentification test for the 2SLS estimator under the assumption of homoskedasticity. L. Hansen (1982) generalized the test to cover the GMM estimator allowing for general heteroskedasticity.

Recall, overidentified models $(\ell>k)$ are special in the sense that there may not be a parameter value $\beta$ such that the moment condition $\mathbb{H}_{0}: \mathbb{E}[Z e]=0$ holds. Thus the model-the overidentifying restrictions - are testable.

For example, take the linear model $Y=\beta_{1}^{\prime} X_{1}+\beta_{2}^{\prime} X_{2}+e$ with $\mathbb{E}\left[X_{1} e\right]=0$ and $\mathbb{E}\left[X_{2} e\right]=0$. It is possible that $\beta_{2}=0$ so that the linear equation may be written as $Y=\beta_{1}^{\prime} X_{1}+e$. However, it is possible that $\beta_{2} \neq 0$. In this case it is impossible to find a value of $\beta_{1}$ such that both $\mathbb{E}\left[X_{1}\left(Y-X_{1}^{\prime} \beta_{1}\right)\right]=0$ and $\mathbb{E}\left[X_{2}\left(Y-X_{1}^{\prime} \beta_{1}\right)\right]=0$ hold simultaneously. In this sense an exclusion restriction can be seen as an overidentifying restriction.

Note that $\bar{g}_{n} \underset{p}{\mathbb{E}}[Z e]$ and thus $\bar{g}_{n}$ can be used to assess the hypothesis $\mathbb{E}[Z e]=0$. Assuming that an efficient weight matrix estimator is used the criterion function at the parameter estimator is $J=$ $J\left(\widehat{\beta}_{\mathrm{gmm}}\right)=n \bar{g}_{n}^{\prime} \widehat{\Omega}^{-1} \bar{g}_{n}$. This is a quadratic form in $\bar{g}_{n}$ and is thus a natural test statistic for $\mathbb{H}_{0}: \mathbb{E}[Z e]=0$. Note that we assume that the criterion function is constructed with an efficient weight matrix estimator. This is important for the distribution theory.

Theorem 13.14 Under Assumption $12.2$ then as $n \rightarrow \infty, J=J\left(\widehat{\beta}_{\mathrm{gmm}}\right) \underset{d}{\rightarrow} \chi_{\ell-k}^{2} \cdot$ For $c$ satisfying $\alpha=1-G_{\ell-k}(c), \mathbb{P}\left[J>c \mid \mathbb{M}_{0}\right] \longrightarrow \alpha$ so the test "Reject $\mathbb{H}_{0}$ if $J>c$ " has asymptotic size $\alpha$. The proof of the theorem is left to Exercise 13.13.

The degrees of freedom of the asymptotic distribution are the number of overidentifying restrictions. If the statistic $J$ exceeds the chi-square critical value we can reject the model. Based on this information alone it is unclear what is wrong but it is typically cause for concern. The GMM overidentification test is a useful by-product of the GMM methodology and it is advisable to report the statistic $J$ whenever GMM is the estimation method. When over-identified models are estimated by GMM it is customary to report the $J$ statistic as a general test of model adequacy.

In Stata the command estat overid afer ivregress gmm can be used to implement the overidentification test. The GMM criterion $J$ and its asymptotic $\mathrm{p}$-value using the $\chi_{\ell-k}^{2}$ distribution are reported.

\subsection{Subset OverIdentification Tests}

In Section $12.32$ we introduced subset overidentification tests for the 2SLS estimator under the assumption of homoskedasticity. In this section we describe how to construct analogous tests for the GMM estimator under general heteroskedasticity.

Recall, subset overidentification tests are used when it is desired to focus attention on a subset of instruments whose validity is questioned. Partition $Z=\left(Z_{a}, Z_{b}\right)$ with dimensions $\ell_{a}$ and $\ell_{b}$, respectively, where $Z_{a}$ contains the instruments which are believed to be uncorrelated with $e$ and $Z_{b}$ contains the instruments which may be correlated with $e$. It is necessary to select this partition so that $\ell_{a}>k$, so that the instruments $Z_{a}$ alone identify the parameters.

Given this partition the maintained hypothesis is $\mathbb{E}\left[Z_{a} e\right]=0$. The null and alternative hypotheses are $\mathbb{H}_{0}: \mathbb{E}\left[Z_{b} e\right]=0$ and $\mathbb{M}_{1}: \mathbb{E}\left[Z_{b} e\right] \neq 0$. The GMM test is constructed as follows. First, estimate the model by efficient GMM with only the smaller set $Z_{a}$ of instruments. Let $\widetilde{J}$ denote the resulting GMM criterion. Second, estimate the model by efficient GMM with the full set $Z=\left(Z_{a}, Z_{b}\right)$ of instruments. Let $\widehat{J}$ denote the resulting GMM criterion. The test statistic is the difference in the criterion functions: $C=\widehat{J}-\widetilde{J}$. This is similar to the GMM distance statistic presented in Section 13.19. The difference is that the distance statistic compares models which differ based on the parameter restrictions while the $C$ statistic compares models based on different instrument sets.

Typically $C \geq 0$. However, this is not necessary and $C<0$ can arise. If this occurs it leads to a nonrejection of $\mathbb{H}_{0}$.

If the smaller instrument set $Z_{a}$ is just-identified so that $\ell_{a}=k$ then $\widetilde{J}=0$ so $C=\widehat{J}$ is simply the standard overidentification test. This is why we have restricted attention to the case $\ell_{a}>k$.

The test has the following large sample distribution.

Theorem 13.15 Under Assumption $12.2$ and $\mathbb{E}\left[Z_{a} X^{\prime}\right]$ has full rank $k$, then as $n \rightarrow \infty, C \rightarrow \underset{d}{\rightarrow} \chi_{\ell_{b}}^{2} .$ For $c$ satisfying $\alpha=1-G_{\ell_{b}}(c), \mathbb{P}\left[C>c \mid \mathbb{H}_{0}\right] \longrightarrow \alpha .$ The test "Reject $\mathbb{H}_{0}$ if $C>c$ " has asymptotic size $\alpha$.

The proof of Theorem $13.15$ is presented in Section 13.28.

In Stata the command estat overid zb afer ivregress gmm can be used to implement a subset overidentification test where $\mathrm{zb}$ is the name(s) of the instruments(s) tested for validity. The statistic $C$ and its asymptotic $p$-value using the $\chi_{\ell_{2}}^{2}$ distribution are reported. 

\subsection{Endogeneity Test}

In Section $12.29$ we introduced tests for endogeneity in the context of 2SLS estimation. Endogeneity tests are simple to implement in the GMM framework as a subset overidentification test. The model is $Y=Z_{1}^{\prime} \beta_{1}+Y_{2}^{\prime} \beta_{2}+e$ where the maintained assumption is that the regressors $Z_{1}$ and excluded instruments $Z_{2}$ are exogenous so that $\mathbb{E}\left[Z_{1} e\right]=0$ and $\mathbb{E}\left[Z_{2} e\right]=0$. The question is whether or not $Y_{2}$ is endogenous. The null hypothesis is $\mathbb{M}_{0}: \mathbb{E}\left[Y_{2} e\right]=0$ with the alternative $\mathbb{H}_{1}: \mathbb{E}\left[Y_{2} e\right] \neq 0$.

The GMM test is constructed as follows. First, estimate the model by efficient GMM using $\left(Z_{1}, Z_{2}\right)$ as instruments for $\left(Z_{1}, Y_{2}\right)$. Let $\widetilde{J}$ denote the resulting GMM criterion. Second, estimate the model by efficient $\mathrm{GMM}^{2}$ using $\left(Z_{1}, Z_{2}, Y_{2}\right)$ as instruments for $\left(Z_{1}, Y_{2}\right)$. Let $\widehat{J}$ denote the resulting GMM criterion. The test statistic is the difference in the criterion functions: $C=\widehat{J}-\widetilde{J}$.

The distribution theory for the test is a special case of overidentification testing.

Theorem $13.16$ Under Assumption $12.2$ and $\mathbb{E}\left[Z_{2} Y_{2}^{\prime}\right]$ has full rank $k_{2}$, then as $n \rightarrow \infty, C \underset{d}{\rightarrow} \chi_{k_{2}}^{2}$. For $c$ satisfying $\alpha=1-G_{k_{2}}(c), \mathbb{P}\left[C>c \mid \mathbb{H}_{0}\right] \rightarrow \alpha$. The test "Reject $\mathbb{H}_{0}$ if $C>c$ " has asymptotic size $\alpha$.

In Stata the command estat endogenous afer ivregress gmm can be used to implement the test for endogeneity. The statistic $C$ and its asymptotic $p$-value using the $\chi_{k_{2}}^{2}$ distribution are reported.

\subsection{Subset Endogeneity Test}

In Section $12.30$ we introduced subset endogeneity tests for 2SLS estimation. GMM tests are simple to implement as subset overidentification tests. The model is $Y=Z_{1}^{\prime} \beta_{1}+Y_{2}^{\prime} \beta_{2}+Y_{3}^{\prime} \beta_{3}+e$ with $\mathbb{E}[Z e]=0$ where the instrument vector is $Z=\left(Z_{1}, Z_{2}\right)$. The $k_{3} \times 1$ variables $Y_{3}$ are treated as endogenous and the $k_{2} \times 1$ variables $Y_{2}$ are treated as potentially endogenous. The hypothesis to test is that $Y_{2}$ is exogenous, or $\mathbb{M}_{0}: \mathbb{E}\left[Y_{2} e\right]=0$ against $\mathbb{H}_{1}: \mathbb{E}\left[Y_{2} e\right] \neq 0$. The test requires that $\ell_{2} \geq\left(k_{2}+k_{3}\right)$ so that the model can be estimated under $\mathbb{H}_{1}$.

The GMM test is constructed as follows. First, estimate the model by efficient GMM using $\left(Z_{1}, Z_{2}\right.$ ) as instruments for $\left(Z_{1}, Y_{2}, Y_{3}\right)$. Let $\widetilde{J}$ denote the resulting GMM criterion. Second, estimate the model by efficient GMM using $\left(Z_{1}, Z_{2}, Y_{2}\right)$ as instruments for $\left(Z_{1}, Y_{2}, Y_{3}\right)$. Let $\widehat{J}$ denote the resulting GMM criterion. The test statistic is the difference in the criterion functions: $C=\widehat{J}-\widetilde{J}$.

The distribution theory for the test is a special case of the theory of overidentification testing.

Theorem $13.17$ Under Assumption $12.2$ and $\mathbb{E}\left[Z_{2}\left(Y_{2}^{\prime}, Y_{3}^{\prime}\right)\right]$ has full rank $k_{2}+k_{3}$, then as $n \rightarrow \infty, C \underset{d}{\longrightarrow} \chi_{k_{2}}^{2}$. For $c$ satisfying $\alpha=1-G_{k_{2}}(c), \mathbb{P}\left[C>c \mid \mathbb{H}_{0}\right] \longrightarrow \alpha$. The test "Reject $\mathbb{H}_{0}$ if $C>c$ " has asymptotic size $\alpha$.

In Stata, the command estat endogenous $\mathrm{x} 2$ afer ivregress gmm can be used to implement the test for endogeneity where $\mathrm{x} 2$ is the name(s) of the variable(s) tested for endogeneity. The statistic $C$ and its asymptotic $\mathrm{p}$-value using the $\chi_{k_{2}}^{2}$ distribution are reported.

${ }^{2}$ If the homoskedastic weight matrix is used this GMM estimator equals least squares, but when the weight matrix allows for heteroskedasticity the efficient GMM estimator does not equal least squares as the model is overidentified. 

\subsection{Nonlinear GMM}

GMM applies whenever an economic or statistical model implies the $\ell \times 1$ moment condition

$$
\mathbb{E}\left[g_{i}(\beta)\right]=0 .
$$

where $g_{i}(\beta)$ is a possibly nonlinear function of the parameters $\beta$. Often, this is all that is known. Identification requires $\ell \geq k=\operatorname{dim}(\beta)$. The GMM estimator minimizes

$$
J(\beta)=n \bar{g}_{n}(\beta)^{\prime} \widehat{\boldsymbol{W}} \bar{g}_{n}(\beta)
$$

for some weight matrix $\widehat{\boldsymbol{W}}$ where

$$
\bar{g}_{n}(\beta)=\frac{1}{n} \sum_{i=1}^{n} g_{i}(\beta) .
$$

The efficient GMM estimator can be constructed by setting

$$
\widehat{\boldsymbol{W}}=\left(\frac{1}{n} \sum_{i=1}^{n} \widehat{g}_{i} \widehat{g}_{i}^{\prime}-\bar{g}_{n} \bar{g}_{n}^{\prime}\right)^{-1},
$$

with $\widehat{g}_{i}=g_{i}(\widetilde{\beta})$ constructed using a preliminary consistent estimator $\widetilde{\beta}$, perhaps obtained with $\widehat{\boldsymbol{W}}=\boldsymbol{I}_{\ell}$. As in the case of the linear model the weight matrix can be iterated until convergence to obtain the iterated GMM estimator.

\section{Proposition 13.1 Distribution of Nonlinear GMM Estimator}

Under general regularity conditions, $\sqrt{n}\left(\widehat{\beta}_{\mathrm{gmm}}-\beta\right) \underset{d}{\longrightarrow} \mathrm{N}\left(0, \boldsymbol{V}_{\beta}\right)$ where

$$
\boldsymbol{V}_{\beta}=\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1}\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \Omega \boldsymbol{W} \boldsymbol{Q}\right)\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1}
$$

with $\Omega=\mathbb{E}\left[g_{i} g_{i}^{\prime}\right]$ and

$$
\boldsymbol{Q}=\mathbb{E}\left[\frac{\partial}{\partial \beta^{\prime}} g_{i}(\beta)\right]
$$

If the efficient weight matrix is used then $\boldsymbol{V}_{\beta}=\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \boldsymbol{Q}\right)^{-1}$.

The proof of this result is omitted as it uses more advanced techniques.

The asymptotic covariance matrices can be estimated by sample counterparts of the population matrices. For the case of a general weight matrix,

$$
\widehat{\boldsymbol{V}}_{\beta}=\left(\widehat{\boldsymbol{Q}}^{\prime} \widehat{\boldsymbol{W}} \widehat{\boldsymbol{Q}}\right)^{-1}\left(\widehat{\boldsymbol{Q}}^{\prime} \widehat{\boldsymbol{W}} \widehat{\Omega} \widehat{\boldsymbol{W}} \widehat{\boldsymbol{Q}}\right)\left(\widehat{\boldsymbol{Q}}^{\prime} \widehat{\boldsymbol{W}} \widehat{\boldsymbol{Q}}\right)^{-1}
$$

where

$$
\begin{gathered}
\widehat{\Omega}=\frac{1}{n} \sum_{i=1}^{n}\left(g_{i}(\widehat{\beta})-\bar{g}\right)\left(g_{i}(\widehat{\beta})-\bar{g}\right)^{\prime} \\
\bar{g}=n^{-1} \sum_{i=1}^{n} g_{i}(\widehat{\beta})
\end{gathered}
$$

and

$$
\widehat{\boldsymbol{Q}}=\frac{1}{n} \sum_{i=1}^{n} \frac{\partial}{\partial \beta^{\prime}} g_{i}(\widehat{\beta}) .
$$

For the case of the iterated efficient weight matrix,

$$
\widehat{\boldsymbol{V}}_{\beta}=\left(\widehat{\boldsymbol{Q}}^{\prime} \widehat{\Omega}^{-1} \widehat{\boldsymbol{Q}}\right)^{-1} .
$$

All of the methods discussed in this chapter - Wald tests, constrained estimation, distance tests, overidentification tests, endogeneity tests - apply similarly to the nonlinear GMM estimator.

\subsection{Bootstrap for GMM}

The bootstrap for 2SLS (Section 12.23) can be used for GMM estimation. The standard bootstrap algorithm generates bootstrap samples by sampling the triplets $\left(Y_{i}^{*}, X_{i}^{*}, Z_{i}^{*}\right)$ independently and with replacement from the original sample. The GMM estimator is applied to the bootstrap sample to obtain the bootstrap estimates $\widehat{\beta}_{\mathrm{gmm}}^{*}$. This is repeated $B$ times to create a sample of $B$ bootstrap draws. Given these draws, bootstrap confidence intervals, including percentile, $\mathrm{BC}$ percentile, $\mathrm{BC}_{a}$ and percentile-t, are calculated conventionally.

For variance and standard error estimation the same cautions apply as for 2SLS. It is difficult to know if the GMM estimator has a finite variance in a given application. It is best to avoid using the bootstrap to calculate standard errors. Instead, use the bootstrap for percentile and percentile-t confidence intervals.

When the model is overidentified, as discussed for 2SLS, bootstrap GMM inference will not achieve an asymptotic refinement unless the bootstrap estimator is recentered to satisfy the orthogonality condition. We now describe the recentering recommended by Hall and Horowitz (1996).

For linear GMM wth weight matrix $\boldsymbol{W}$ the recentered GMM bootstrap estimator is

$$
\widehat{\beta}_{\mathrm{gmm}}^{* *}=\left(\boldsymbol{X}^{* \prime} \boldsymbol{Z}^{*} \boldsymbol{W}^{*} \boldsymbol{Z}^{* \prime} \boldsymbol{X}^{*}\right)^{-1}\left(\boldsymbol{X}^{* \prime} \boldsymbol{Z}^{*} \boldsymbol{W}^{*}\left(\boldsymbol{Z}^{* \prime} \boldsymbol{Y}^{*}-\boldsymbol{Z}^{\prime} \widehat{\boldsymbol{e}}\right)\right)
$$

where $\boldsymbol{W}^{*}$ is the bootstrap version of $\boldsymbol{W}$ and $\widehat{\boldsymbol{e}}=\boldsymbol{Y}-\boldsymbol{X} \widehat{\beta}_{\mathrm{gmm}}$. For efficient GMM,

$$
\boldsymbol{W}^{*}=\left(\frac{1}{n} \sum_{i=1}^{n} Z_{i}^{*} Z_{i}^{* \prime}\left(Y_{i}^{*}-X_{i}^{* \prime} \widetilde{\beta}^{*}\right)^{2}\right)^{-1}
$$

for preliminary estimator $\widetilde{\beta}^{*}$.

For nonlinear GMM (Section 13.25) the bootstrap criterion function is modified. The recentered bootstrap criterion is

$$
\begin{aligned}
&J^{* *}(\beta)=n\left(\bar{g}_{n}^{*}(\beta)-\bar{g}_{n}\left(\widehat{\beta}_{\mathrm{gmm}}\right)\right)^{\prime} \boldsymbol{W}^{*}\left(\bar{g}_{n}^{*}(\beta)-\bar{g}_{n}\left(\widehat{\beta}_{\mathrm{gmm}}\right)\right) \\
&\bar{g}_{n}^{*}(\beta)=\frac{1}{n} \sum_{i=1}^{n} g_{i}^{*}(\beta)
\end{aligned}
$$

where $\bar{g}_{n}\left(\widehat{\beta}_{\mathrm{gmm}}\right)$ is from the sample not from the bootstrap data. The bootstrap estimator is

$$
\widehat{\beta}_{\mathrm{gmm}}^{* *}=\operatorname{argmin} J^{* *}(\beta) .
$$

The bootstrap can be used to calculate the p-value of the GMM overidentification test. For the GMM estimator with an efficient weight matrix the standard overidentification test is the Hansen $J$ statistic

$$
J=n \bar{g}_{n}\left(\widehat{\beta}_{\mathrm{gmm}}\right)^{\prime} \widehat{\Omega}^{-1} \bar{g}_{n}\left(\widehat{\beta}_{\mathrm{gmm}}\right) .
$$

The recentered bootstrap analog is

$$
J^{* *}=n\left(\bar{g}_{n}^{*}\left(\widehat{\beta}_{\mathrm{gmm}}^{* *}\right)-\bar{g}_{n}\left(\widehat{\beta}_{\mathrm{gmm}}\right)\right)^{\prime} \widehat{\Omega}^{*-1}\left(\bar{g}_{n}^{*}\left(\widehat{\beta}_{\mathrm{gmm}}^{* *}\right)-\bar{g}_{n}\left(\widehat{\beta}_{\mathrm{gmm}}\right)\right)
$$

On each bootstrap sample $J^{* *}(b)$ is calculated and stored. The bootstrap p-value is

$$
p^{*}=\frac{1}{B} \sum_{b=1}^{B} \mathbb{1}\left\{J^{* *}(b)>S\right\} .
$$

This bootstrap p-value is asymptotically valid since $J^{* *}$ satisfies the overidentified moment conditions.

\subsection{Conditional Moment Equation Models}

In many contexts, an economic model implies conditional moment restriction of the form

$$
\mathbb{E}\left[e_{i}(\beta) \mid Z_{i}\right]=0
$$

where $e_{i}(\beta)$ is some $s \times 1$ function of the observation and the parameters. In many cases $s=1$. It turns out that this conditional moment restriction is more powerful than the unconditional moment equation model discussed throughout this chapter.

For example, the linear model $Y=X^{\prime} \beta+e$ with instruments $Z$ falls into this class under the assumption $\mathbb{E}[e \mid Z]=0$. In this case $e_{i}(\beta)=Y_{i}-X_{i}^{\prime} \beta$.

It is also helpful to realize that conventional regression models also fall into this class except that in this case $X=Z$. For example, in linear regression $e_{i}(\beta)=Y_{i}-X_{i}^{\prime} \beta$, while in a nonlinear regression model $e_{i}(\beta)=Y_{i}-m\left(X_{i}, \beta\right)$. In a joint model of the conditional expectation $\mathbb{E}[Y \mid X=x]=x^{\prime} \beta$ and variance $\operatorname{var}[Y \mid X=x]=f(x)^{\prime} \gamma$, then

$$
e_{i}(\beta, \gamma)=\left\{\begin{array}{c}
Y_{i}-X_{i}^{\prime} \beta \\
\left(Y_{i}-X_{i}^{\prime} \beta\right)^{2}-f\left(X_{i}\right)^{\prime} \gamma
\end{array} .\right.
$$

Here $s=2$.

Given a conditional moment restriction an unconditional moment restriction can always be constructed. That is for any $\ell \times 1$ function $\phi(Z, \beta)$ we can set $g_{i}(\beta)=\phi\left(Z_{i}, \beta\right) e_{i}(\beta)$ which satisfies $\mathbb{E}\left[g_{i}(\beta)\right]=$ 0 and hence defines an unconditional moment equation model. The obvious problem is that the class of functions $\phi$ is infinite. Which should be selected?

This is equivalent to the problem of selection of the best instruments. If $Z \in \mathbb{R}$ is a valid instrument satisfying $\mathbb{E}[e \mid Z]=0$, then $Z, Z^{2}, Z^{3}$,..., etc., are all valid instruments. Which should be used?

One solution is to construct an infinite list of potent instruments and then use the first $\ell$. How is $\ell$ to be determined? This is an area of theory still under development. One study of this problem is Donald and Newey (2001).

Another approach is to construct the optimal instrument which minimizes the asymptotic variance. The form was uncovered by Chamberlain (1987). Take the case $s=1$. Let

$$
R_{i}=\mathbb{E}\left[\frac{\partial}{\partial \beta} e_{i}(\beta) \mid Z_{i}\right]
$$

and $\sigma_{i}^{2}=\mathbb{E}\left[e_{i}(\beta)^{2} \mid Z_{i}\right]$. Then the optimal instrument is $A_{i}=-\sigma_{i}^{-2} R_{i}$. The optimal moment is $g_{i}(\beta)=$ $A_{i} e_{i}(\beta)$. Setting $g_{i}(\beta)$ to be this choice (which is $k \times 1$, so is just-identified) yields the GMM estimator with lowest asymptotic variance. In practice $A_{i}$ is unknown, but its form helps us think about construction of good instruments. In the linear model $e_{i}(\beta)=Y_{i}-X_{i}^{\prime} \beta$ note that $R_{i}=-\mathbb{E}\left[X_{i} \mid Z_{i}\right]$ and $\sigma_{i}^{2}=\mathbb{E}\left[e_{i}^{2} \mid Z_{i}\right]$. This means the optimal instrument is $A_{i}=\sigma_{i}^{-2} \mathbb{E}\left[X_{i} \mid Z_{i}\right]$. In the case of linear regression $X_{i}=Z_{i}$ so $A_{i}=\sigma_{i}^{-2} Z_{i}$. Hence efficient GMM is equivalent to GLS!

In the case of endogenous variables note that the efficient instrument $A_{i}$ involves the estimation of the conditional mean of $X$ given $Z$. In other words, to get the best instrument for $X$ we need the best conditional mean model for $X$ given $Z$ not just an arbitrary linear projection. The efficient instrument is also inversely proportional to the conditional variance of $e$. This is the same as the GLS estimator; namely that improved efficiency can be obtained if the observations are weighted inversely to the conditional variance of the errors.

\subsection{Technical Proofs*}

Proof of Theorem 13.12 Set $\widetilde{e}_{i}=Y_{i}-X_{i}^{\prime} \widehat{\beta}_{\mathrm{cgmm}}$ and $\widehat{e}_{i}=Y_{i}-X_{i}^{\prime} \widehat{\beta}_{\mathrm{gmm}}$. By standard covariance matrix analysis $\widehat{\Omega} \vec{p} \Omega$ and $\widetilde{\Omega} \vec{p} \Omega$. Thus we can replace $\widehat{\Omega}$ and $\widetilde{\Omega}$ in the criteria without affecting the asymptotic distribution. In particular

$$
\begin{aligned}
\widetilde{J}\left(\widehat{\beta}_{\mathrm{cgmm}}\right) &=\frac{1}{n} \widetilde{\boldsymbol{e}}^{\prime} \boldsymbol{Z} \widetilde{\Omega}^{-1} \boldsymbol{Z}^{\prime} \widetilde{\boldsymbol{e}} \\
&=\frac{1}{n} \widetilde{\boldsymbol{e}}^{\prime} \boldsymbol{Z} \widehat{\Omega}^{-1} \boldsymbol{Z}^{\prime} \widetilde{\boldsymbol{e}}+o_{p}(1) .
\end{aligned}
$$

Now observe that

$$
\boldsymbol{Z}^{\prime} \widetilde{\boldsymbol{e}}=\boldsymbol{Z}^{\prime} \widehat{\boldsymbol{e}}-\boldsymbol{Z}^{\prime} \boldsymbol{X}\left(\widehat{\beta}_{\mathrm{cgmm}}-\widehat{\beta}_{\mathrm{gmm}}\right) .
$$

Thus

$$
\begin{aligned}
\frac{1}{n} \widetilde{\boldsymbol{e}}^{\prime} \boldsymbol{Z} \widehat{\Omega}^{-1} \boldsymbol{Z}^{\prime} \widetilde{\boldsymbol{e}} &=\frac{1}{n} \widehat{\boldsymbol{e}}^{\prime} \boldsymbol{Z} \widehat{\Omega}^{-1} \boldsymbol{Z}^{\prime} \widehat{\boldsymbol{e}}-\frac{2}{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\widehat{\beta}_{\mathrm{gmm}}\right)^{\prime} \boldsymbol{X}^{\prime} \boldsymbol{Z} \widehat{\Omega}^{-1} \boldsymbol{Z}^{\prime} \widehat{\boldsymbol{e}} \\
&+\frac{1}{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\widehat{\beta}_{\mathrm{gmm}}\right)^{\prime} \boldsymbol{X}^{\prime} \boldsymbol{Z} \widehat{\Omega}^{-1} \boldsymbol{Z}^{\prime} \boldsymbol{X}\left(\widehat{\beta}_{\mathrm{cgmm}}-\widehat{\beta}_{\mathrm{gmm}}\right) \\
&=\widehat{J}\left(\widehat{\beta}_{\mathrm{gmm}}\right)+\frac{1}{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\widehat{\beta}_{\mathrm{gmm}}\right)^{\prime} \boldsymbol{X}^{\prime} \boldsymbol{Z} \widehat{\Omega}^{-1} \boldsymbol{Z}^{\prime} \boldsymbol{X}\left(\widehat{\beta}_{\mathrm{cgmm}}-\widehat{\beta}_{\mathrm{gmm}}\right)
\end{aligned}
$$

where the second equality holds because $\boldsymbol{X}^{\prime} \boldsymbol{Z} \widehat{\Omega}^{-1} Z^{\prime} \widehat{\boldsymbol{e}}=0$ is the first-order condition for $\widehat{\beta}_{\mathrm{gmm}}$. By (13.16) and Theorem 13.4, under $\mathbb{M}_{0}$

$$
\begin{aligned}
\sqrt{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\widehat{\beta}_{\mathrm{gmm}}\right) &=-\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \Omega^{-1} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{R}\left(\boldsymbol{R}^{\prime}\left(\boldsymbol{X}^{\prime} \boldsymbol{Z} \Omega^{-1} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime} \sqrt{n}\left(\widehat{\beta}_{\mathrm{gmm}}-\beta\right)+o_{p}(1) \\
& \underset{d}{\longrightarrow}\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \boldsymbol{Q}\right)^{-1} \boldsymbol{R} Z
\end{aligned}
$$

where

$$
\begin{aligned}
Z & \sim \mathrm{N}\left(0, \boldsymbol{V}_{\boldsymbol{R}}\right) \\
\boldsymbol{V}_{\boldsymbol{R}} &=\left(\boldsymbol{R} \boldsymbol{V}^{\prime}\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \boldsymbol{Q}\right)^{-1} \boldsymbol{R}\right)^{-1} .
\end{aligned}
$$

Putting together (13.25), (13.26), (13.27) and (13.28),

$$
\begin{aligned}
D &=\widetilde{J}\left(\widehat{\beta}_{\mathrm{cgmm}}\right)-\widehat{J}\left(\widehat{\beta}_{\mathrm{gmm}}\right) \\
&=\sqrt{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\widehat{\beta}_{\mathrm{gmm}}\right)^{\prime} \frac{1}{n} \boldsymbol{X}^{\prime} \boldsymbol{Z} \widehat{\Omega}^{-1} \frac{1}{n} \boldsymbol{Z}^{\prime} \boldsymbol{X} \sqrt{n}\left(\widehat{\beta}_{\mathrm{cgmm}}-\widehat{\beta}_{\mathrm{gmm}}\right) \\
& \underset{d}{\longrightarrow} Z^{\prime} \boldsymbol{V}_{\boldsymbol{R}}^{-1} Z \sim \chi_{q}^{2}
\end{aligned}
$$

because $V_{R}>0$ and $\mathrm{Z}$ is $q \times 1$.

Proof of Theorem 13.15 Let $\widetilde{\beta}$ denote the GMM estimator obtained with the instrument set $Z_{a}$ and let $\widehat{\beta}$ denote the GMM estimator obtained with the instrument set $Z$. Set $\widetilde{e}_{i}=Y_{i}-X_{i}^{\prime} \widetilde{\beta}, \widehat{e}{ }_{i}=Y_{i}-X_{i}^{\prime} \widehat{\beta}$,

$$
\begin{aligned}
&\widetilde{\Omega}=n^{-1} \sum_{i=1}^{n} Z_{a i} Z_{a i}^{\prime} \widetilde{e}_{i}^{2} \\
&\widehat{\Omega}=n^{-1} \sum_{i=1}^{n} Z_{i} Z_{i}^{\prime} \widehat{e}_{i}^{2}
\end{aligned}
$$

Let $\boldsymbol{R}$ be the $\ell \times \ell_{a}$ selector matrix so that $Z_{a}=\boldsymbol{R}^{\prime} Z$. Note that

$$
\widetilde{\Omega}=\boldsymbol{R}^{\prime} n^{-1} \sum_{i=1}^{n} Z_{i} Z_{i}^{\prime} \widetilde{e}_{i}^{2} \boldsymbol{R} .
$$

By standard covariance matrix analysis, $\widehat{\Omega} \underset{p}{\rightarrow} \Omega$ and $\widetilde{\Omega} \underset{p}{\rightarrow} \boldsymbol{R}^{\prime} \Omega \boldsymbol{R}$. Also, $\frac{1}{n} \boldsymbol{Z}^{\prime} \boldsymbol{X} \underset{p}{\rightarrow} \boldsymbol{Q}$, say. By the CLT, $n^{-1 / 2} \boldsymbol{Z}^{\prime} \boldsymbol{e} \underset{d}{\longrightarrow} Z$ where $Z \sim \mathrm{N}(0, \Omega)$. Then

$$
\begin{aligned}
n^{-1 / 2} \boldsymbol{Z}^{\prime} \widehat{\boldsymbol{e}} &=\left(\boldsymbol{I}_{\ell}-\left(\frac{1}{n} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)\left(\frac{1}{n} \boldsymbol{X}^{\prime} \boldsymbol{Z} \widehat{\Omega}^{-1} \frac{1}{n} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1}\left(\frac{1}{n} \boldsymbol{X}^{\prime} \boldsymbol{Z}\right) \widehat{\Omega}^{-1}\right) n^{-1 / 2} \boldsymbol{Z}^{\prime} \boldsymbol{e} \\
& \rightarrow\left(\boldsymbol{I}_{\ell}-\boldsymbol{Q}\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \boldsymbol{Q}\right)^{-1} \boldsymbol{Q}^{\prime} \Omega^{-1}\right) Z
\end{aligned}
$$

and

$$
\begin{aligned}
n^{-1 / 2} \boldsymbol{Z}_{a}^{\prime} \widetilde{\boldsymbol{e}} &=\boldsymbol{R}^{\prime}\left(\boldsymbol{I}_{\ell}-\left(\frac{1}{n} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)\left(\frac{1}{n} \boldsymbol{X}^{\prime} \boldsymbol{Z} \boldsymbol{R}^{-1} \boldsymbol{R}^{\prime} \frac{1}{n} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)^{-1}\left(\frac{1}{n} \boldsymbol{X}^{\prime} \boldsymbol{Z}\right) \boldsymbol{R}^{-1} \boldsymbol{R}^{\prime}\right) n^{-1 / 2} \boldsymbol{Z}^{\prime} \boldsymbol{e} \\
& \underset{d}{\longrightarrow} \boldsymbol{R}^{\prime}\left(\boldsymbol{I}_{\ell}-\boldsymbol{Q}\left(\boldsymbol{Q}^{\prime} \boldsymbol{R}\left(\boldsymbol{R}^{\prime} \Omega \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime} \boldsymbol{Q}\right)^{-1} \boldsymbol{Q}^{\prime} \boldsymbol{R}\left(\boldsymbol{R}^{\prime} \Omega \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime}\right) Z
\end{aligned}
$$

jointly.

By linear rotations of $Z$ and $\boldsymbol{R}$ we can set $\Omega=\boldsymbol{I}_{\ell}$ to simplify the notation. Thus setting $\boldsymbol{P}_{\boldsymbol{Q}}=\boldsymbol{Q}\left(\boldsymbol{Q}^{\prime} \boldsymbol{Q}\right)^{-1} \boldsymbol{Q}^{\prime}$, $\boldsymbol{P}_{\boldsymbol{R}}=\boldsymbol{R}\left(\boldsymbol{R}^{\prime} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime}$ and $Z \sim \mathrm{N}\left(0, \boldsymbol{I}_{\ell}\right)$ we have

$$
\widehat{J} \underset{d}{\longrightarrow} Z^{\prime}\left(I_{\ell}-\boldsymbol{P}_{\mathbf{Q}}\right) Z
$$

and

$$
\widetilde{J} \underset{d}{\rightarrow} Z^{\prime}\left(\boldsymbol{P}_{\boldsymbol{R}}-\boldsymbol{P}_{\boldsymbol{R}} \boldsymbol{Q}\left(\boldsymbol{Q}^{\prime} \boldsymbol{P}_{R} \boldsymbol{Q}\right)^{-1} \boldsymbol{Q}^{\prime} \boldsymbol{P}_{\boldsymbol{R}}\right) Z
$$

It follows that

$$
C=\widehat{J}-\widetilde{J} \underset{d}{\longrightarrow} \mathrm{Z}^{\prime} A \mathrm{Z}
$$

where

$$
\boldsymbol{A}=\left(\boldsymbol{I}_{\ell}-\boldsymbol{P}_{Q}-\boldsymbol{P}_{\boldsymbol{R}}+\boldsymbol{P}_{R} \boldsymbol{Q}\left(\boldsymbol{Q}^{\prime} \boldsymbol{P}_{R} \boldsymbol{Q}\right)^{-1} \boldsymbol{Q}^{\prime} \boldsymbol{P}_{R}\right) .
$$

This is a quadratic form in a standard normal vector and the matrix $\boldsymbol{A}$ is idempotent (this is straightforward to check). $Z^{\prime} A Z$ is thus distributed $\chi_{d}^{2}$ with degrees of freedom $d$ equal to

$$
\begin{aligned}
\operatorname{rank}(\boldsymbol{A}) &=\operatorname{tr}\left(\boldsymbol{I}_{\ell}-\boldsymbol{P}_{\boldsymbol{Q}}-\boldsymbol{P}_{\boldsymbol{R}}+\boldsymbol{P}_{\boldsymbol{R}} \boldsymbol{Q}\left(\boldsymbol{Q}^{\prime} \boldsymbol{P}_{\boldsymbol{R}} \boldsymbol{Q}\right)^{-1} \boldsymbol{Q}^{\prime} \boldsymbol{P}_{\boldsymbol{R}}\right) \\
&=\ell-k-\ell_{a}+k=\ell_{b}
\end{aligned}
$$

Thus the asymptotic distribution of $C$ is $\chi_{\ell_{b}}^{2}$ as claimed. 

\subsection{Exercises}

Exercise 13.1 Take the model

$$
\begin{aligned}
Y &=X^{\prime} \beta+e \\
\mathbb{E}[X e] &=0 \\
e^{2} &=Z^{\prime} \gamma+\eta \\
\mathbb{E}[Z \eta] &=0 .
\end{aligned}
$$

Find the method of moments estimators $(\widehat{\beta}, \widehat{\gamma})$ for $(\beta, \gamma)$

Exercise 13.2 Take the model $Y=X^{\prime} \beta+e$ with $\mathbb{E}[e \mid Z]=0$. Let $\widehat{\beta}_{\mathrm{gmm}}$ be the GMM estimator using the weight matrix $\boldsymbol{W}_{n}=\left(\boldsymbol{Z}^{\prime} \boldsymbol{Z}\right)^{-1}$. Under the assumption $\mathbb{E}\left[e^{2} \mid Z\right]=\sigma^{2}$ show that

$$
\sqrt{n}(\widehat{\beta}-\beta) \underset{d}{\longrightarrow} \mathrm{N}\left(0, \sigma^{2}\left(\boldsymbol{Q}^{\prime} \boldsymbol{M}^{-1} \boldsymbol{Q}\right)^{-1}\right)
$$

where $\boldsymbol{Q}=\mathbb{E}\left[Z X^{\prime}\right]$ and $\boldsymbol{M}=\mathbb{E}\left[Z Z^{\prime}\right]$

Exercise 13.3 Take the model $Y=X^{\prime} \beta+e$ with $\mathbb{E}[Z e]=0$. Let $\widetilde{e}_{i}=Y_{i}-X_{i}^{\prime} \widetilde{\beta}$ where $\widetilde{\beta}$ is consistent for $\beta$ (e.g. a GMM estimator with some weight matrix). An estimator of the optimal GMM weight matrix is

$$
\widehat{\boldsymbol{W}}=\left(\frac{1}{n} \sum_{i=1}^{n} Z_{i} Z_{i}^{\prime} \widetilde{e}_{i}^{2}\right)^{-1} .
$$

Show that $\widehat{\boldsymbol{W}} \underset{p}{\longrightarrow} \Omega^{-1}$ where $\Omega=\mathbb{E}\left[Z Z^{\prime} e^{2}\right]$.

Exercise $13.4$ In the linear model estimated by GMM with general weight matrix $\boldsymbol{W}$ the asymptotic variance of $\widehat{\beta}_{\mathrm{gmm}}$ is

$$
\boldsymbol{V}=\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1} \boldsymbol{Q}^{\prime} \boldsymbol{W} \Omega \boldsymbol{W} \boldsymbol{Q}\left(\boldsymbol{Q}^{\prime} \boldsymbol{W} \boldsymbol{Q}\right)^{-1}
$$

(a) Let $\boldsymbol{V}_{0}$ be this matrix when $\boldsymbol{W}=\Omega^{-1}$. Show that $\boldsymbol{V}_{0}=\left(\boldsymbol{Q}^{\prime} \Omega^{-1} \boldsymbol{Q}\right)^{-1}$.

(b) We want to show that for any $\boldsymbol{W}, \boldsymbol{V}-\boldsymbol{V}_{0}$ is positive semi-definite (for then $\boldsymbol{V}_{0}$ is the smaller possible covariance matrix and $W=\Omega^{-1}$ is the efficient weight matrix). To do this start by finding matrices $\boldsymbol{A}$ and $\boldsymbol{B}$ such that $\boldsymbol{V}=\boldsymbol{A}^{\prime} \Omega \boldsymbol{A}$ and $\boldsymbol{V}_{0}=\boldsymbol{B}^{\prime} \Omega \boldsymbol{B}$.

(c) Show that $\boldsymbol{B}^{\prime} \Omega \boldsymbol{A}=\boldsymbol{B}^{\prime} \Omega \boldsymbol{B}$ and therefore that $\boldsymbol{B}^{\prime} \Omega(\boldsymbol{A}-\boldsymbol{B})=0$.

(d) Use the expressions $\boldsymbol{V}=\boldsymbol{A}^{\prime} \mathbf{\Omega} \boldsymbol{A}, \boldsymbol{A}=\boldsymbol{B}+(\boldsymbol{A}-\boldsymbol{B})$, and $\boldsymbol{B}^{\prime} \boldsymbol{\Omega}(\boldsymbol{A}-\boldsymbol{B})=0$ to show that $\boldsymbol{V} \geq \boldsymbol{V}_{0}$.

Exercise $13.5$ Prove Theorem 13.8.

Exercise 13.6 Derive the constrained GMM estimator (13.16).

Exercise 13.7 Show that the constrained GMM estimator (13.16) with the efficient weight matrix is (13.19).

Exercise $13.8$ Prove Theorem 13.9.

Exercise 13.9 Prove Theorem 13.10. Exercise $13.10$ The equation of interest is $Y=m(X, \beta)+e$ with $\mathbb{E}[Z e]=0$ where $m(x, \beta)$ is a known function, $\beta$ is $k \times 1$ and $Z$ is $\ell \times 1$. Show how to construct an efficient GMM estimator for $\beta$.

Exercise 13.11 As a continuation of Exercise $12.7$ derive the efficient GMM estimator using the instrument $Z=\left(\begin{array}{ll}X & X^{2}\end{array}\right)^{\prime}$. Does this differ from 2SLS and/or OLS?

Exercise 13.12 In the linear model $Y=X^{\prime} \beta+e$ with $\mathbb{E}[X e]=0$ the GMM criterion function for $\beta$ is

$$
J(\beta)=\frac{1}{n}(\boldsymbol{Y}-\boldsymbol{X} \beta)^{\prime} \boldsymbol{X} \widehat{\Omega}^{-1} \boldsymbol{X}^{\prime}(\boldsymbol{Y}-\boldsymbol{X} \beta)
$$

where $\widehat{\Omega}=n^{-1} \sum_{i=1}^{n} X_{i} X_{i}^{\prime} \widehat{e}_{i}^{2}, \widehat{e}_{i}=Y_{i}-X_{i}^{\prime} \widehat{\beta}$ are the OLS residuals, and $\widehat{\beta}=\left(\boldsymbol{X}^{\prime} \boldsymbol{X}\right)^{-1} \boldsymbol{X}^{\prime} \boldsymbol{Y}$ is least squares. The GMM estimator of $\beta$ subject to the restriction $r(\beta)=0$ is

$$
\widetilde{\beta}=\underset{r(\beta)=0}{\operatorname{argmin}} J_{n}(\beta) .
$$

The GMM test statistic (the distance statistic) of the hypothesis $r(\beta)=0$ is

$$
D=J(\tilde{\beta})=\min _{r(\beta)=0} J(\beta) .
$$

(a) Show that you can rewrite $J(\beta)$ in (13.29) as

$$
J(\beta)=n(\beta-\widehat{\beta})^{\prime} \widehat{\boldsymbol{V}}_{\beta}^{-1}(\beta-\widehat{\beta})
$$

and thus $\widetilde{\beta}$ is the same as the minimum distance estimator.

(b) Show that under linear hypotheses the distance statistic $D$ in (13.30) equals the Wald statistic.

Exercise 13.13 Take the linear model $Y=X^{\prime} \beta+e$ with $\mathbb{E}[Z e]=0$. Consider the GMM estimator $\widehat{\beta}$ of $\beta$. Let $J=n \bar{g}_{n}(\widehat{\beta})^{\prime} \widehat{\Omega}^{-1} \bar{g}_{n}(\widehat{\beta})$ denote the test of overidentifying restrictions. Show that $J \underset{d}{\longrightarrow} \chi_{\ell-k}^{2}$ as $n \rightarrow \infty$ by demonstrating each of the following.

(a) Since $\Omega>0$, we can write $\Omega^{-1}=\boldsymbol{C} \boldsymbol{C}^{\prime}$ and $\Omega=\boldsymbol{C}^{\prime-1} \boldsymbol{C}^{-1}$ for some matrix $\boldsymbol{C}$.

(b) $J=n\left(\boldsymbol{C}^{\prime} \bar{g}_{n}(\widehat{\beta})\right)^{\prime}\left(\boldsymbol{C}^{\prime} \widehat{\Omega} \boldsymbol{C}\right)^{-1} \boldsymbol{C}^{\prime} \bar{g}_{n}(\widehat{\beta})$.

(c) $\boldsymbol{C}^{\prime} \bar{g}_{n}(\widehat{\beta})=\boldsymbol{D}_{n} \boldsymbol{C}^{\prime} \bar{g}_{n}(\beta)$ where $\bar{g}_{n}(\beta)=\frac{1}{n} \boldsymbol{Z}^{\prime} \boldsymbol{e}$ and

$$
\boldsymbol{D}_{n}=\boldsymbol{I}_{\ell}-\boldsymbol{C}^{\prime}\left(\frac{1}{n} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)\left(\left(\frac{1}{n} \boldsymbol{X}^{\prime} \boldsymbol{Z}\right) \widehat{\Omega}^{-1}\left(\frac{1}{n} \boldsymbol{Z}^{\prime} \boldsymbol{X}\right)\right)^{-1}\left(\frac{1}{n} \boldsymbol{X}^{\prime} \boldsymbol{Z}\right) \widehat{\Omega}^{-1} \boldsymbol{C}^{\prime-1}
$$

(d) $\boldsymbol{D}_{n} \underset{p}{\longrightarrow} \boldsymbol{I}_{\ell}-\boldsymbol{R}\left(\boldsymbol{R}^{\prime} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime}$ where $\boldsymbol{R}=\boldsymbol{C}^{\prime} \mathbb{E}\left[Z X^{\prime}\right]$.

(e) $n^{1 / 2} \boldsymbol{C}^{\prime} \bar{g}_{n}(\beta) \underset{d}{\longrightarrow} u \sim \mathrm{N}\left(0, \boldsymbol{I}_{\ell}\right)$.

(f) $J \underset{d}{\longrightarrow} u^{\prime}\left(I_{\ell}-\boldsymbol{R}\left(\boldsymbol{R}^{\prime} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime}\right) u$.

(g) $u^{\prime}\left(\boldsymbol{I}_{\ell}-\boldsymbol{R}\left(\boldsymbol{R}^{\prime} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime}\right) u \sim \chi_{\ell-k}^{2}$.

Hint: $\boldsymbol{I}_{\ell}-\boldsymbol{R}\left(\boldsymbol{R}^{\prime} \boldsymbol{R}\right)^{-1} \boldsymbol{R}^{\prime}$ is a projection matrix. Exercise 13.14 Take the model $Y=X^{\prime} \beta+e$ with $\mathbb{E}[Z e]=0, Y \in \mathbb{R}, X \in \mathbb{R}^{k}, Z \in \mathbb{R}^{\ell}, \ell \geq k$. Consider the statistic

$$
\begin{aligned}
J(\beta) &=n \bar{m}_{n}(\beta)^{\prime} \boldsymbol{W} \bar{m}_{n}(\beta) \\
\bar{m}_{n}(\beta) &=\frac{1}{n} \sum_{i=1}^{n} Z_{i}\left(Y_{i}-X_{i}^{\prime} \beta\right)
\end{aligned}
$$

for some weight matrix $W>0$.

(a) Take the hypothesis $\mathbb{I}_{0}: \beta=\beta_{0}$. Derive the asymptotic distribution of $J\left(\beta_{0}\right)$ under $\mathbb{H}_{0}$ as $n \rightarrow \infty$.

(b) What choice for $W$ yields a known asymptotic distribution in part (a)? (Be specific about degrees of freedom.)

(c) Write down an appropriate estimator $\widehat{\boldsymbol{W}}$ for $W$ which takes advantage of $\mathbb{M}_{0}$. (You do not need to demonstrate consistency or unbiasedness.)

(d) Describe an asymptotic test of $\mathbb{H}_{0}$ against $\mathbb{M}_{1}: \beta \neq \beta_{0}$ based on this statistic.

(e) Use the result in part (d) to construct a confidence region for $\beta$. What can you say about the form of this region? For example, does the confidence region take the form of an ellipse, similar to conventional confidence regions?

Exercise 13.15 Consider the model $Y=X^{\prime} \beta+e$ with $\mathbb{E}[Z e]=0$ and

$$
\boldsymbol{R}^{\prime} \beta=0
$$

with $Y \in \mathbb{R}, X \in \mathbb{R}^{k}, Z \in \mathbb{R}^{\ell}, \ell>k$. The matrix $\boldsymbol{R}$ is $k \times q$ with $1 \leq q<k$. You have a random sample $\left(Y_{i}, X_{i}, Z_{i}: i=1, \ldots, n\right)$.

For simplicity, assume the efficient weight matrix $\boldsymbol{W}=\left(\mathbb{E}\left[Z Z^{\prime} e^{2}\right]\right)^{-1}$ is known.

(a) Write out the GMM estimator $\widehat{\beta}$ ignoring constraint (13.31).

(b) Write out the GMM estimator $\widetilde{\beta}$ adding the constraint (13.31).

(c) Find the asymptotic distribution of $\sqrt{n}(\widetilde{\beta}-\beta)$ as $n \rightarrow \infty$ under Assumption (13.31).

Exercise $13.16$ The observed data is $\left\{Y_{i}, X_{i}, Z_{i}\right\} \in \mathbb{R} \times \mathbb{R}^{k} \times \mathbb{R}^{\ell}, k>1$ and $\ell>k>1, i=1, \ldots, n$. The model is $Y=X^{\prime} \beta+e$ with $\mathbb{E}[Z e]=0$.

(a) Given a weight matrix $\boldsymbol{W}>0$ write down the GMM estimator $\widehat{\beta}$ for $\beta$.

(b) Suppose the model is misspecified. Specifically, assume that for some $\delta \neq 0$,

$$
\begin{aligned}
e &=\delta n^{-1 / 2}+u \\
\mathbb{E}[u \mid Z] &=0
\end{aligned}
$$

with $\mu_{Z}=\mathbb{E}[Z] \neq 0$. Show that (13.32) implies that $\mathbb{E}[Z e] \neq 0$.

(c) Express $\sqrt{n}(\widehat{\beta}-\beta)$ as a function of $\boldsymbol{W}, n, \delta$, and the variables $\left(X_{i}, Z_{i}, u_{i}\right)$.

(d) Find the asymptotic distribution of $\sqrt{n}(\widehat{\beta}-\beta)$ under Assumption (13.32). Exercise $13.17$ The model is $Y=Z \beta+X \gamma+e$ with $\mathbb{E}[e \mid Z]=0, X \in \mathbb{R}$ and $Z \in \mathbb{R}$. $X$ is potentially endogenous and $Z$ is exogenous. Someone suggests estimating $(\beta, \gamma)$ by GMM using the pair $\left(Z, Z^{2}\right)$ as instruments. Is this feasible? Under what conditions is this a valid estimator?

Exercise $13.18$ The observations are i.i.d., $\left(Y_{i}, X_{i}, Q_{i}: i=1, \ldots, n\right)$, where $X$ is $k \times 1$ and $Q$ is $m \times 1$. The model is $Y=X^{\prime} \beta+e$ with $\mathbb{E}[X e]=0$ and $\mathbb{E}[Q e]=0$. Find the efficient GMM estimator for $\beta$.

Exercise 13.19 You want to estimate $\mu=\mathbb{E}[Y]$ under the assumption that $\mathbb{E}[X]=0$, where $Y$ and $X$ are scalar and observed from a random sample. Find an efficient GMM estimator for $\mu$.

Exercise 13.20 Consider the model $Y=X^{\prime} \beta+e$ given $\mathbb{E}[Z e]=0$ and $\boldsymbol{R}^{\prime} \beta=0$. The dimensions are $X \in R^{k}$ and $Z \in R^{\ell}$ with $\ell>k$. The matrix $\boldsymbol{R}$ is $k \times q, 1 \leq q<k$. Derive an efficient GMM estimator for $\beta$.

Exercise 13.21 Take the linear equation $Y=X^{\prime} \beta+e$ and consider the following estimators of $\beta$.

1. $\widehat{\beta}$ : 2SLS using the instruments $Z_{1}$.

2. $\widetilde{\beta}: 2$ SLS using the instruments $Z_{2}$.

3. $\bar{\beta}$ : GMM using the instruments $Z=\left(Z_{1}, Z_{2}\right)$ and the weight matrix

$$
\boldsymbol{W}=\left(\begin{array}{cc}
\left(\boldsymbol{Z}_{1}^{\prime} \boldsymbol{Z}_{1}\right)^{-1} \lambda & 0 \\
0 & \left(\boldsymbol{Z}_{2}^{\prime} \boldsymbol{Z}_{2}\right)^{-1}(1-\lambda)
\end{array}\right)
$$

for $\lambda \in(0,1)$.

Find an expression for $\bar{\beta}$ which shows that it is a specific weighted average of $\widehat{\beta}$ and $\widetilde{\beta}$.

Exercise 13.22 Consider the just-identified model $Y=X_{1}^{\prime} \beta_{1}+X_{2}^{\prime} \beta_{2}+e$ with $\mathbb{E}[Z e]=0$ where $X=\left(X_{1}^{\prime}\right.$ $\left.X_{2}^{\prime}\right)^{\prime} \in \mathbb{R}^{k}$ and $Z \in \mathbb{R}^{k}$. We want to test $\mathbb{H}_{0}: \beta_{1}=0$. Three econometricians are called for advice.

- Econometrician 1 proposes testing $\mathbb{M}_{0}$ by a Wald statistic.

- Econometrician 2 suggests testing $\mathbb{M}_{0}$ by the GMM Distance Statistic.

- Econometrician 3 suggests testing $\mathbb{M}_{0}$ using the test of overidentifying restrictions.

You are asked to settle this dispute. Explain the advantages and/or disadvantages of the different procedures in this specific context.

Exercise 13.23 Take the model $Y=X^{\prime} \beta+e$ with $\mathbb{E}[X e]=0$ and $\beta=\boldsymbol{Q} \theta$, where $\beta$ is $k \times 1, \boldsymbol{Q}$ is $k \times m$ with $m<k, \boldsymbol{Q}$ is known, and $\theta$ is $m \times 1$. The observations $\left(Y_{i}, X_{i}\right)$ are i.i.d. across $i=1, \ldots, n$.

Under these assumptions what is the efficient estimator of $\theta$ ?

Exercise 13.24 Take the model $Y=\theta+e$ with $\mathbb{E}[X e]=0, Y \in \mathbb{R}, X \in \mathbb{R}^{k}$ and $\left(Y_{i}, X_{i}\right)$ a random sample.

(a) Find the efficient GMM estimator of $\theta$.

(b) Is this model over-identified or just-identified?

(c) Find the GMM test statistic for over-identification. Exercise 13.25 Take the model $Y=X^{\prime} \beta+e$ with $\mathbb{E}[X e]=0$ where $X$ contains an intercept so $\mathbb{E}[e]=0$. An enterprising econometrician notices that this implies the $n$ moment conditions

$$
\mathbb{E}\left[e_{i}\right]=0, i=1, \ldots, n .
$$

Given an $n \times n$ weight matrix $\boldsymbol{W}$, this implies a GMM criterion

$$
J(\beta)=(\boldsymbol{Y}-\boldsymbol{X} \beta)^{\prime} \boldsymbol{W}(\boldsymbol{Y}-\boldsymbol{X} \beta) .
$$

(a) Under i.i.d. sampling, show that the efficient weight matrix is $\boldsymbol{W}=\sigma^{-2} \boldsymbol{I}_{n}$ where $\sigma^{2}=\mathbb{E}\left[e^{2}\right]$.

(b) Using the weight matrix $\boldsymbol{W}=\sigma^{-2} \boldsymbol{I}_{n}$ find the GMM estimator $\widehat{\beta}$ that minimizes $J(\beta)$.

(c) Find a simple expression for the minimized criteria $J(\widehat{\beta})$.

(d) Theorem $13.14$ says that criterion such as $J(\widehat{\beta})$ are asymptotically $\chi_{\ell-k}^{2}$ where $\ell$ is the number of moments. While the assumptions of Theorem $13.14$ do not apply to this context, what is $\ell$ here? That is, which $\chi^{2}$ distribution is the asserted asymptotic distribution?

(e) Does the answer in (d) make sense? Explain your reasoning.

Exercise 13.26 Take the model $Y=X^{\prime} \beta+e$ with $\mathbb{E}[e \mid X]=0$ and $\mathbb{E}\left[e^{2} \mid X\right]=\sigma^{2}$. An econometrician more enterprising than the one in previous question notices that this implies the $n k$ moment conditions

$$
\mathbb{E}\left[X_{i} e_{i}\right]=0, i=1, \ldots, n .
$$

We can write the moments using matrix notation as $\mathbb{E}\left[\bar{X}^{\prime}(\boldsymbol{Y}-\boldsymbol{X} \beta)\right]$ where

$$
\overline{\boldsymbol{X}}=\left(\begin{array}{cccc}
X_{1}^{\prime} & 0 & \cdots & 0 \\
0 & X_{2}^{\prime} & & 0 \\
\vdots & \vdots & & \vdots \\
0 & 0 & \cdots & X_{n}^{\prime}
\end{array}\right) \text {. }
$$

Given an $n k \times n k$ weight matrix $\boldsymbol{W}$ this implies a GMM criterion

$$
J(\beta)=(\boldsymbol{Y}-\boldsymbol{X} \beta)^{\prime} \overline{\boldsymbol{X}} \boldsymbol{W} \overline{\boldsymbol{X}}^{\prime}(\boldsymbol{Y}-\boldsymbol{X} \beta) .
$$

(a) Calculate $\Omega=\mathbb{E}\left[\overline{\boldsymbol{X}}^{\prime} \boldsymbol{e} \boldsymbol{e}^{\prime} \overline{\boldsymbol{X}}\right]$.

(b) The econometrician decides to set $\boldsymbol{W}=\Omega^{-}$, the Moore-Penrose generalized inverse of $\Omega$. (See Section A.6.) Note: A useful fact is that for a vector $\boldsymbol{a},\left(\boldsymbol{a} \boldsymbol{a}^{\prime}\right)^{-}=\boldsymbol{a} \boldsymbol{a}^{\prime}\left(\boldsymbol{a}^{\prime} \boldsymbol{a}\right)^{-2}$.

(c) Find the GMM estimator $\widehat{\beta}$ that minimizes $J(\beta)$.

(d) Find a simple expression for the minimized criterion $J(\widehat{\beta})$.

(e) Comment on whether the $\chi^{2}$ approximation from Theorem $13.14$ is appropriate for $J(\widehat{\beta})$.

Exercise 13.27 Continuation of Exercise 12.22, based on the empirical work reported in Acemoglu, Johnson, and Robinson (2001).

(a) Re-estimate the model estimated in part (j) by efficient GMM. Use the 2SLS estimates as the firststep for the weight matrix and then calculate the GMM estimator using this weight matrix without further iteration. Report the estimates and standard errors. (b) Calculate and report the $J$ statistic for overidentification.

(c) Compare the GMM and 2SLS estimates. Discuss your findings.

Exercise 13.28 Continuation of Exercise 12.24, which involved estimation of a wage equation by 2 SLS.

(a) Re-estimate the model in part (a) by efficient GMM. Do the results change meaningfully?

(b) Re-estimate the model in part (d) by efficient GMM. Do the results change meaningfully?

(c) Report the $J$ statistic for overidentification.