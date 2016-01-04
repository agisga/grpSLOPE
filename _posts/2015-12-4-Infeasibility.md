---
layout: post
title: Derivation of infeasibility as stopping criterion
---

As one stopping criterion for the Group SLOPE optimization problem, [Brzyski et. al. (2015) *Group SLOPE &mdash; adaptive selection of groups of predictors*](http://arxiv.org/abs/1511.09078) derives in Appendix I a measure of infeasibility of the following form:

$$
\begin{equation}
\label{infeas}
\mathrm{infeas}(\mu^k) = \max \left\\{ J^D\subscript{\lambda, I}\left( \tilde{\xi} \right) - 1, 0\right\\},
\end{equation}
$$

where $J^D\subscript{\lambda, I}$ denotes the dual norm to the Group Sorted L1 norm $J\subscript{\lambda, I}$.

That is, the infeasibility is measured based on the unit ball of the dual norm to the Group Sorted L1 norm, and it is greater than zero if and only if $\tilde{\xi}$ is outside that unit ball. However, it is not clear how the dual norm in question can be easily evaluated. In the following, we will construct an infeasibility criterion which is easy to compute.

Let $w\subscript{\mathbb{I}} := (\\|w\subscript{I\subscript{1}}\\|\subscript{2}, \ldots, \\|w\subscript{I\subscript{J}}\\|\subscript{2})^T$ denotes the vector of Euclidean norms of the groups of elements in $w$, which are given by the partition $I$, and let $J\subscript{\lambda}$ denote the Sorted L1 norm. Then the dual norm of the Group Sorted L1 norm is given by,

$$J^D\subscript{\lambda, I}(x) = \max\subscript{w}\left\\{x^T w : J\subscript{\lambda}(w\subscript{\mathbb{I}}) = 1\right\\},$$

Using the Cauchy-Schwarz inequality we have that

$$x^T w = \sum\subscript{i} w\subscript{I\subscript{i}}^T x\subscript{I\subscript{i}}
\leq \sum\subscript{i} \left\\|w\subscript{I\subscript{i}}\right\\|\subscript{2} \left\\|x\subscript{I\subscript{i}}\right\\|\subscript{2}.$$

It follows that

$$J^D\subscript{\lambda, I}(x) \leq \max\subscript{u}\left\\{\sum u\subscript{i} \\|x\subscript{I\subscript{i}}\\|\subscript{2} : J\subscript{\lambda}(u) = 1\right\\}.$$

In fact, equality holds in the above expression. Assume that $\tilde{u}$ is the maximum:

$$\tilde{u} = \mathrm{argmax}\subscript{u}\left\\{\sum u\subscript{i} \\|x\subscript{I\subscript{i}}\\|\subscript{2} : J\subscript{\lambda}(u) = 1\right\\}.$$

Define $\tilde{w}$ by

$$\tilde{w}\subscript{I\subscript{i}} := \tilde{u}\subscript{i} \frac{x\subscript{I\subscript{i}}}{\left\\|x\subscript{I\subscript{i}}\right\\|\subscript{2}}.$$

Then it holds that

$$
\begin{eqnarray}
x^T \tilde{w} &=& \sum\subscript{i} \tilde{w}\subscript{I\subscript{i}}^T x\subscript{I\subscript{i}} \nonumber \\\\
&=& \sum\subscript{i} \frac{\tilde{u}\subscript{i}}{\left\\|x\subscript{I\subscript{i}}\right\\|\subscript{2}} x\subscript{I\subscript{i}}^T x\subscript{I\subscript{i}} \nonumber \\\\
&=& \sum\subscript{i} \tilde{u}\subscript{i} \left\\|x\subscript{I\subscript{i}}\right\\|\subscript{2}, \nonumber \\\\
\end{eqnarray}
$$

and

$$J\subscript{\lambda}(\tilde{w}\subscript{\mathbb{I}}) = \sum \lambda\subscript{i} \left\lVert \tilde{w}\subscript{I\subscript{(i)}} \right\rVert\subscript{2}
= \sum \lambda\subscript{i} \left\lvert \tilde{u} \right\rvert\subscript{(i)} = J\subscript{\lambda}(\tilde{u}) = 1.$$

Thus, we conclude that

$$J^D\subscript{\lambda, I}(x) = \max\subscript{u}\left\\{\sum u\subscript{i} \\|x\subscript{I\subscript{i}}\\|\subscript{2} : J\subscript{\lambda}(u) = 1\right\\} = J^D\subscript{\lambda}(x\subscript{\mathbb{I}}),$$

where $x\subscript{\mathbb{I}} := (\\|x\subscript{I\subscript{1}}\\|\subscript{2}, \ldots, \\|x\subscript{I\subscript{J}}\\|\subscript{2})^T$ denotes the vector of Euclidean norms of the groups of elements in $x$, and where $J^D\subscript{\lambda}$ denotes the dual norm to the Sorted L1 norm.

Therefore, the unit ball of the norm $J^D\subscript{\lambda, I}(x)$ is the set $\left\\{x\in\mathbb{R}^n : J^D\subscript{\lambda}(x\subscript{\mathbb{I}}) \leq 1 \right\\}$, i.e., the set of all $x$ such that $x\subscript{\mathbb{I}}$ is in the unit ball of the dual norm to the Sorted L1 norm.

Proposition 1.1 in [Bogdan et. al. (2013), *Statistical estimation and testing via the sorted L1 norm*](http://arxiv.org/abs/1310.1969) establishes that a vector $u$ is contained in the unit ball of the dual norm to the Sorted L1 norm, if and only if

$$\forall j : \sum\subscript{i=1}^j \left\lvert x \right\rvert\subscript{(i)} \leq \sum\subscript{i=1}^j \lambda\subscript{i}.$$

Based on $(\ref{infeas})$ and the above considerations about the unit ball of the dual norm, a measure of infeasibility for the Group SLOPE problem is given by

$$\max \left\\{ 0, \max\subscript{j}\left\\{\sum\subscript{i=1}^j \left( \left\lVert x\subscript{I\subscript{(i)}} \right\rVert\subscript{2} - \lambda\subscript{i} \right) \right\\} \right\\}.$$
