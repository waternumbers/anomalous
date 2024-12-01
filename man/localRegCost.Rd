% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loc_reg.R
\name{localRegCost}
\alias{localRegCost}
\title{R6 class for local regression Cost Functions}
\description{
VERY experimental - do not use
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-localRegCost-length}{\code{localRegCost$length()}}
\item \href{#method-localRegCost-new}{\code{localRegCost$new()}}
\item \href{#method-localRegCost-baseCost}{\code{localRegCost$baseCost()}}
\item \href{#method-localRegCost-pointCost}{\code{localRegCost$pointCost()}}
\item \href{#method-localRegCost-collectiveCost}{\code{localRegCost$collectiveCost()}}
\item \href{#method-localRegCost-param}{\code{localRegCost$param()}}
\item \href{#method-localRegCost-clone}{\code{localRegCost$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-localRegCost-length"></a>}}
\if{latex}{\out{\hypertarget{method-localRegCost-length}{}}}
\subsection{Method \code{length()}}{
Get the length of time series
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{localRegCost$length()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-localRegCost-new"></a>}}
\if{latex}{\out{\hypertarget{method-localRegCost-new}{}}}
\subsection{Method \code{new()}}{
Initialise the cost function
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{localRegCost$new(x, family = c("gaussian", "symmetric"))}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{observations as for gauss_reg}

\item{\code{family}}{for fitting}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-localRegCost-baseCost"></a>}}
\if{latex}{\out{\hypertarget{method-localRegCost-baseCost}{}}}
\subsection{Method \code{baseCost()}}{
Compute the non-anomalous cost of a segment
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{localRegCost$baseCost(a, b, pen = 0)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{a}}{start of period}

\item{\code{b}}{end of period}

\item{\code{pen}}{penalty cost}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-localRegCost-pointCost"></a>}}
\if{latex}{\out{\hypertarget{method-localRegCost-pointCost}{}}}
\subsection{Method \code{pointCost()}}{
Compute the point anomaly cost of a time step
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{localRegCost$pointCost(b, pen)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{b}}{time step}

\item{\code{pen}}{penalty cost}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-localRegCost-collectiveCost"></a>}}
\if{latex}{\out{\hypertarget{method-localRegCost-collectiveCost}{}}}
\subsection{Method \code{collectiveCost()}}{
Compute the non-anomalous cost of a segment
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{localRegCost$collectiveCost(a, b, pen, len)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{a}}{start of period}

\item{\code{b}}{end of period}

\item{\code{pen}}{penalty cost}

\item{\code{len}}{minimum number of observations}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-localRegCost-param"></a>}}
\if{latex}{\out{\hypertarget{method-localRegCost-param}{}}}
\subsection{Method \code{param()}}{
Compute parameters of a segment if anomalous
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{localRegCost$param(a, b)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{a}}{start of period}

\item{\code{b}}{end of period}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-localRegCost-clone"></a>}}
\if{latex}{\out{\hypertarget{method-localRegCost-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{localRegCost$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
