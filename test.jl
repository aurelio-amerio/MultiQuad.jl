using FastGaussQuadrature
using Unitful
using Base.Threads
nthreads()
nodes, weights = gausslegendre( 10)

length(nodes)


a = 1u"km"
unit(a)