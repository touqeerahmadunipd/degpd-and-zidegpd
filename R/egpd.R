## extended generalised Pareto negative log-likelihood functions

## model 1 ##

# F.expr <- expression((1 - (1 + xi * x / exp(lpsi))^(-1/xi))^exp(lkappa))
# D(F.expr, "x")
# f.expr <- expression(-((1 - (1 + xi * x/exp(lpsi))^(-1/xi))^(exp(lkappa) - 1) * (exp(lkappa) * 
#     ((1 + xi * x/exp(lpsi))^((-1/xi) - 1) * ((-1/xi) * (xi/exp(lpsi)))))))
# nll.expr <- expression(-log(-((1 - (1 + xi * x/exp(lpsi))^(-1/xi))^(exp(lkappa) - 1) * (exp(lkappa) * 
#     ((1 + xi * x/exp(lpsi))^((-1/xi) - 1) * ((-1/xi) * (xi/exp(lpsi))))))))

# egpd1d0 <- function(pars, X1, X2, X3, y, dupid, duplicate) {
#   # vectors of parameters
#   lpsi <- X1 %*% pars[[1]]
#   xi <- X2 %*% pars[[2]]
#   lkappa <- X3 %*% pars[[3]]
#   # interim stuff
#   .e1 <- 1/xi
#   .e4 <- y * xi/exp(lpsi)
#   if (any(.e4 <= -1)) return(1e20)
#   # neg log lik vector
#   out <- -((exp(lkappa) - 1) * log(1 - 1/(1 + .e4)^.e1) + lkappa - ((1 + .e1) * log1p(.e4) + lpsi))
#   # overall neg log lik
#   sum(out)
# }
# 
# egpd1d12 <- function(pars, X1, X2, X3, y, dupid, duplicate) {
#   # vectors of parameters
#   lpsi <- X1 %*% pars[[1]]
#   xi <- X2 %*% pars[[2]]
#   lkappa <- X3 %*% pars[[3]]
#   # interim stuff
#   .e1 <- exp(lpsi)
#   .e2 <- y * xi
#   .e3 <- .e2/.e1
#   .e4 <- 1/xi
#   .e5 <- 1 + .e3
#   .e6 <- 1 + .e4
#   .e7 <- .e5^.e4
#   .e8 <- exp(lkappa)
#   .e9 <- 1 - 1/.e7
#   .e10 <- .e5^.e6
#   .e11 <- .e8 - 1
#   .e12 <- .e9^.e11
#   .e13 <- log1p(.e3)
#   .e14 <- .e8 - 2
#   .e15 <- .e4 + 2
#   .e16 <- .e9^.e14
#   .e17 <- .e5^.e15
#   .e18 <- .e17 * .e1
#   .e19 <- xi^2
#   .e20 <- 2 * .e6
#   .e21 <- y * .e6
#   .e22 <- log(.e9)
#   .e23 <- .e10 * .e1
#   .e24 <- .e5^.e20
#   .e25 <- xi * .e7
#   .e26 <- .e16 * .e11
#   .e27 <- xi * .e10
#   .e28 <- .e5^(.e4 - 1)
#   .e30 <- .e13/.e25 - y/.e23
#   .e31 <- .e2 * .e6
#   .e32 <- .e19 * .e10
#   .e34 <- .e7 * .e13/xi
#   .e35 <- .e24 * .e1
#   .e36 <- 1/.e10
#   .e38 <- y * .e28/.e1
#   .e41 <- .e36 - .e31/.e18
#   .e42 <- 2 * .e11
#   .e44 <- .e13/.e32 - .e21/.e18
#   .e47 <- .e21 * .e7/.e1 - .e10 * .e13/.e19
#   .e48 <- .e38 - .e34
#   .e49 <- .e12 * .e41
#   .e51 <- .e12 * .e8 * .e22
#   .e52 <- .e12/.e10
#   .e53 <- .e9^(.e8 - (2 + .e42))
#   .e54 <- .e49 + y * .e16 * .e11/.e35
#   .e57 <- .e51/.e10 + .e52
#   .e60 <- .e26 * .e30/.e27 - .e12 * .e44
#   .e61 <- 2/xi
#   .e62 <- .e5^.e61
#   .e63 <- .e18^2
#   .e64 <- .e26 * .e22
#   .e65 <- .e9^(.e8 - 3)
#   .e66 <- .e9^(.e8 - (1 + .e42))
#   .e67 <- .e65 * .e14
#   .e68 <- .e5^(1 - .e4)
#   .e69 <- .e5^(.e4 - .e20)
#   .e70 <- .e54 * .e53
#   .e71 <- .e57 * .e53
#   .e72 <- .e60 * .e53
#   .e73 <- .e23^2
#   .e74 <- .e35^2
#   .e75 <- .e12 * .e22
#   .e76 <- .e16 * .e41
#   .e77 <- .e16 * .e44
#   .e78 <- .e16 + .e64
#   .e79 <- .e67 * .e30
#   .e80 <- .e18 - .e2 * .e10 * .e15
#   .e83 <- .e5^(.e20 - 1)
#   .e84 <- .e27^2
#   .e85 <- .e25^2
#   .e86 <- .e32^2
#   .e90 <- y * .e10 * .e15/.e1 - .e17 * .e13/.e19
#   .e93 <- xi * .e12 * .e6 * .e69
#   .e94 <- xi * .e62
#   .e95 <- xi * .e47
#   # final matrix
#   out <- matrix(0, length(y), 9)
#   # first derivatives
#   out[, 1] <- .e54 * .e10/.e12 # d1
#   out[, 2] <- .e60 * .e10/.e12 # d2
#   out[, 3] <- -(.e57 * .e10/.e12) # d3
#   # second derivatives
#   out[, 4] <- y * (.e70 * .e11/.e1 + (.e10 * (xi * (.e80/.e63 + .e69/.e1) * 
#             .e12 * .e6 - ((.e76/.e10 + y * .e65 * .e14/(.e5^(.e6 + 
#             .e20) * .e1))/.e1 + (.e35 - 2 * (.e31 * .e83)) * 
#             .e16/.e74) * .e11) - xi * .e54 * .e6 * .e7/.e1)/.e12) # d11
#   out[, 5] <- y * (((((((.e23 - .e31 * .e7)/.e73 + (xi * .e28 * 
#             .e13/.e85 - .e36)/.e1) * .e16 - .e79/.e23)/.e27 + 
#             (.e77/.e10 + .e19 * .e16 * .e6 * .e7 * .e30/.e84)/.e1) * 
#             .e11 - (.e80 * .e6/.e63 + (xi^3 * .e6 * .e7 * .e13/.e86 - 
#             1/(xi * .e17))/.e1) * .e12) * .e10 - xi * .e60 * 
#             .e6 * .e7/.e1)/.e12 + .e72 * .e11/.e1) # d12
#   out[, 6] <- -(y * 
#             (.e71 * .e11 + (.e10 * (.e8 * (.e93 * .e22 - (.e64/.e10 + 
#                 .e16/.e10)/.e10) + .e93 - .e26/.e24) - xi * .e57 * 
#                 .e6 * .e7)/.e12)/.e1) # d13
#   out[, 7] <- (((((.e16 * (y * (1/(.e27 * 
#         .e1) + .e1 * .e47/.e73) - (.e7 + .e38 - .e34) * .e13/.e85) + 
#         .e79 * .e48/.e94)/.e10 - .e77 * .e48/.e62)/xi - (.e10 + 
#         .e95) * .e16 * .e30/.e84) * .e11 - .e12 * (y * (.e6 * 
#         .e1 * .e90/.e63 + 2/(.e19 * .e17 * .e1)) - xi * (2 * 
#         .e10 + .e95) * .e13/.e86)) * .e10 + .e60 * .e47)/.e12 - 
#         .e72 * .e68 * .e11 * .e48/xi # d22
#   out[, 8] <- -(((((.e64/.e62 + 
#         .e16/.e62) * .e48/.e27 - .e75 * .e47/.e24) * .e8 + .e26 * 
#         .e48/(xi * .e5^(1 + 3/xi)) - .e12 * .e47/.e24) * .e10 + 
#         .e57 * .e47)/.e12 - .e71 * .e68 * .e11 * .e48/xi) # d23
#   out[, 9] <- -((((.e12 + 
#         .e51)/.e10 + .e52) * .e10/.e12 - .e57 * .e66 * .e10) * 
#         .e8 * .e22) # d33
#   out
# }
# 
# egpd1d34 <- function(pars, X1, X2, X3, y, dupid, duplicate) {
#   # vectors of parameters
#   lpsi <- X1 %*% pars[[1]]
#   xi <- X2 %*% pars[[2]]
#   lkappa <- X3 %*% pars[[3]]
#   # interim stuff
#     .e1 <- exp(lpsi)
#   .e2 <- y * xi
#   .e3 <- .e2/.e1
#   .e4 <- 1/xi
#   .e5 <- 1 + .e3
#   .e6 <- 1 + .e4
#   .e7 <- .e5^.e4
#   .e8 <- exp(lkappa)
#   .e9 <- 1 - 1/.e7
#   .e10 <- .e5^.e6
#   .e11 <- log1p(.e3)
#   .e12 <- .e8 - 1
#   .e13 <- .e8 - 2
#   .e14 <- .e4 + 2
#   .e15 <- .e9^.e13
#   .e16 <- xi^2
#   .e17 <- .e5^.e14
#   .e18 <- 2 * .e6
#   .e19 <- .e9^.e12
#   .e20 <- .e4 - 1
#   .e21 <- .e17 * .e1
#   .e22 <- .e5^.e20
#   .e23 <- y * .e6
#   .e24 <- log(.e9)
#   .e25 <- .e10 * .e1
#   .e27 <- y * .e22/.e1
#   .e29 <- .e7 * .e11/xi
#   .e30 <- .e5^.e18
#   .e31 <- 2/xi
#   .e32 <- xi * .e7
#   .e33 <- .e27 - .e29
#   .e34 <- .e10 * .e11
#   .e36 <- .e23 * .e7/.e1
#   .e38 <- .e36 - .e34/.e16
#   .e39 <- .e5^.e31
#   .e40 <- 2 * .e12
#   .e41 <- .e15 * .e12
#   .e42 <- .e2 * .e6
#   .e43 <- xi * .e10
#   .e44 <- .e8 - 3
#   .e45 <- .e16 * .e10
#   .e46 <- .e9^.e44
#   .e48 <- .e11/.e32 - y/.e25
#   .e49 <- 1/.e10
#   .e50 <- .e30 * .e1
#   .e52 <- .e8 - (2 + .e40)
#   .e53 <- .e49 - .e42/.e21
#   .e55 <- .e11/.e45 - .e23/.e21
#   .e56 <- .e21^2
#   .e57 <- .e9^.e52
#   .e58 <- .e41 * .e24
#   .e59 <- .e46 * .e13
#   .e60 <- .e4 - .e18
#   .e61 <- xi * .e39
#   .e62 <- .e5^.e60
#   .e63 <- .e18 - 1
#   .e64 <- .e25^2
#   .e66 <- .e19 * .e8 * .e24
#   .e67 <- .e32^2
#   .e68 <- .e19/.e10
#   .e69 <- .e5^.e63
#   .e70 <- .e19 * .e53
#   .e71 <- y * .e10
#   .e72 <- .e71 * .e14
#   .e73 <- .e17 * .e11
#   .e74 <- .e72/.e1
#   .e75 <- xi * .e38
#   .e76 <- .e73/.e16
#   .e77 <- .e74 - .e76
#   .e78 <- .e21 - .e2 * .e10 * .e14
#   .e79 <- .e45^2
#   .e81 <- .e8 - (1 + .e40)
#   .e82 <- .e43^2
#   .e83 <- 1 - .e4
#   .e84 <- .e70 + y * .e15 * .e12/.e50
#   .e85 <- .e19 * .e55
#   .e87 <- .e66/.e10 + .e68
#   .e89 <- .e41 * .e48/.e43
#   .e90 <- .e89 - .e85
#   .e91 <- .e7 + .e27
#   .e92 <- .e50^2
#   .e93 <- .e9^.e81
#   .e94 <- .e15 + .e58
#   .e95 <- 2 * .e10
#   .e96 <- xi * .e22
#   .e97 <- .e15 * .e53
#   .e98 <- .e5^.e83
#   .e99 <- .e42 * .e7
#   .e100 <- 3/xi
#   .e101 <- .e19 * .e24
#   .e102 <- .e15 * .e55
#   .e103 <- .e59 * .e48
#   .e104 <- .e25 - .e99
#   .e105 <- .e91 - .e29
#   .e106 <- xi * .e19
#   .e107 <- xi * .e6
#   .e108 <- .e5^(.e4 - 2)
#   .e109 <- .e10 + .e75
#   .e110 <- .e30 * .e11
#   .e111 <- 1 + .e100
#   .e112 <- .e18 + .e31
#   .e113 <- 2 * (.e23 * .e69/.e1)
#   .e114 <- .e106 * .e6
#   .e115 <- .e43 * .e1
#   .e116 <- .e22 * .e11
#   .e119 <- .e113 - 2 * (.e110/.e16)
#   .e120 <- .e114 * .e62
#   .e123 <- .e5^.e111
#   .e127 <- y * (1/.e115 + .e1 * .e38/.e64) - .e105 * .e11/.e67
#   .e129 <- .e104/.e64 + (.e96 * .e11/.e67 - .e49)/.e1
#   .e130 <- .e15/.e39
#   .e131 <- .e5^.e112
#   .e132 <- .e6 + .e18
#   .e133 <- .e95 + .e75
#   .e134 <- xi^3
#   .e135 <- xi * .e17
#   .e137 <- .e78/.e56 + .e62/.e1
#   .e139 <- .e58/.e39 + .e130
#   .e140 <- .e15/.e10
#   .e141 <- .e59 * .e33
#   .e142 <- .e5^.e132
#   .e144 <- .e16 * .e17 * .e1
#   .e146 <- .e58/.e10 + .e140
#   .e147 <- .e50 - 2 * (.e42 * .e69)
#   .e148 <- .e38/.e30
#   .e149 <- 1/.e135
#   .e150 <- xi * .e131
#   .e151 <- xi * .e133
#   .e152 <- .e129 * .e15
#   .e155 <- .e78 * .e6/.e56 + (.e134 * .e6 * .e7 * .e11/.e79 - .e149)/.e1
#   .e156 <- .e19 + .e66
#   .e161 <- .e148 + y * (1/.e21 - .e107 * .e1 * .e77/.e56)
#   .e163 <- y * (.e6 * .e1 * .e77/.e56 + 2/.e144) - .e151 * .e11/.e79
#   .e164 <- .e15 * .e127
#   .e165 <- .e142 * .e1
#   .e166 <- 1 + .e18
#   .e167 <- xi * .e123
#   .e168 <- .e150 * .e1
#   .e170 <- .e156/.e10 + .e68
#   .e171 <- .e147 * .e15
#   .e172 <- .e164 + .e103 * .e33/.e61
#   .e173 <- xi * .e137
#   .e174 <- .e152 - .e103/.e25
#   .e175 <- .e9^(.e8 - (.e40 + 3))
#   .e176 <- .e70 * .e24
#   .e177 <- .e101 * .e55
#   .e178 <- .e97/.e10
#   .e179 <- .e15 * .e119
#   .e180 <- .e102/.e10
#   .e181 <- .e59 * .e24
#   .e182 <- .e16 * .e15
#   .e183 <- .e155 * .e19
#   .e184 <- .e146/.e10
#   .e185 <- .e172/.e10
#   .e188 <- .e94 * .e48/.e43 - .e177
#   .e189 <- .e109 * .e15
#   .e190 <- .e161 * .e19
#   .e191 <- .e176 + y * .e94/.e50
#   .e192 <- .e101 * .e38
#   .e194 <- .e97 * .e33/.e61
#   .e195 <- .e179 * .e1
#   .e196 <- .e41 * .e33
#   .e197 <- .e41/.e30
#   .e198 <- .e102 * .e33
#   .e199 <- .e180 + .e182 * .e6 * .e7 * .e48/.e82
#   .e200 <- .e116/.e16
#   .e202 <- y * .e46 * .e13
#   .e205 <- y * .e108 * .e20/.e1
#   .e206 <- .e2 * .e108
#   .e208 <- .e173 * .e19 * .e6
#   .e209 <- .e120 * .e24
#   .e210 <- .e174/.e43
#   .e211 <- (.e185 - .e198/.e39)/xi
#   .e212 <- (.e178 + .e202/.e165)/.e1
#   .e213 <- .e199/.e1
#   .e215 <- .e189 * .e48/.e82
#   .e216 <- (.e116 + .e96)/xi
#   .e217 <- .e171/.e92
#   .e219 <- .e19 * .e38/.e30
#   .e220 <- .e192/.e30
#   .e224 <- .e7 + .e29
#   .e227 <- .e8 * (.e209 - .e184) + .e120 - .e197
#   .e228 <- y * (.e141/.e168 - .e195/.e92)
#   .e231 <- y * (.e205 - .e200)/.e1
#   .e233 <- .e206 * .e20/.e1
#   .e235 <- (.e210 + .e213) * .e12 - .e183
#   .e237 <- (.e211 - .e215) * .e12 - .e19 * .e163
#   .e241 <- (.e139 * .e33/.e43 - .e220) * .e8 + .e196/.e167 - .e219
#   .e242 <- .e84 * .e57
#   .e244 <- (.e194 + .e228) * .e12 - .e190
#   .e245 <- .e90 * .e57
#   .e246 <- (.e6 * .e7 * .e11 + .e7)/xi
#   .e247 <- .e216 - (.e22 + .e233)
#   .e249 <- .e34/xi + .e95
#   .e250 <- .e7/xi
#   .e251 <- .e5^(.e31 - 1)
#   .e254 <- y * .e7/.e1
#   .e255 <- .e231 - (.e11 * (.e27 - .e224)/xi + .e27)/xi
#   .e256 <- .e208 - (.e212 + .e217) * .e12
#   .e257 <- .e9^(.e8 - 4)
#   .e258 <- .e6 * .e33
#   .e259 <- .e87 * .e57
#   .e260 <- .e246 - .e91 * .e6
#   .e262 <- .e257 * .e44
#   .e264 <- .e5^.e166
#   .e265 <- .e5^(1 + .e31)
#   .e267 <- .e1^2
#   .e270 <- y * (.e258 - .e250)/.e1 - (.e11 * (.e36 - .e249/xi) + .e254)/xi
#   .e271 <- .e57 * .e98
#   .e272 <- .e5^(.e4 - .e166)
#   .e273 <- xi * .e87
#   .e274 <- xi * .e15
#   .e275 <- .e87 * .e93
#   .e276 <- .e7 * .e1
#   .e279 <- 2 * (.e39 * .e11/xi)
#   .e280 <- 2 * (y * .e251/.e1)
#   .e281 <- xi * .e84
#   .e282 <- xi * .e90
#   .e283 <- .e274 * .e6
#   .e284 <- .e170 * .e57
#   .e285 <- .e188 * .e57
#   .e286 <- .e129 * .e46
#   .e287 <- .e155 * .e15
#   .e288 <- .e191 * .e57
#   .e290 <- .e84 * .e175 * .e52
#   .e292 <- .e87 * .e175 * .e52
#   .e294 <- .e90 * .e175 * .e52
#   .e295 <- .e161 * .e15
#   .e296 <- .e15 * .e163
#   .e297 <- .e15 * .e24
#   .e299 <- .e46 * .e53 * .e13
#   .e300 <- .e59 * .e55
#   .e302 <- .e181/.e10 + .e46/.e10
#   .e304 <- .e181/.e39 + .e46/.e39
#   .e305 <- .e46 * .e127
#   .e306 <- .e46 + .e181
#   .e307 <- .e262 * .e48
#   .e309 <- .e10 * .e227 - .e273 * .e6 * .e7
#   .e310 <- .e5^(.e18 - 2)
#   .e311 <- .e5^(4 * .e6)
#   .e312 <- .e5^(4/xi)
#   .e313 <- .e61^2
#   .e314 <- 2 * .e7
#   .e315 <- .e99/.e1
#   .e317 <- .e235 * .e10 - .e282 * .e6 * .e7/.e1
#   .e319 <- .e237 * .e10 + .e90 * .e38
#   .e321 <- .e241 * .e10 + .e87 * .e38
#   .e323 <- .e170 * .e10/.e19
#   .e325 <- .e244 * .e10 + .e84 * .e38
#   .e326 <- .e247 * .e15
#   .e327 <- .e84 * .e93
#   .e328 <- .e275 * .e10
#   .e329 <- .e90 * .e93
#   .e336 <- .e271 * .e255 - .e271 * .e33/xi
#   .e337 <- .e15 * (.e280 - .e279)
#   .e339 <- (1 - 2 * (.e78 * .e17 * .e1/.e56)) * .e1 * .e77
#   .e342 <- .e98 * .e11/.e16 + y * .e83/.e276
#   .e343 <- .e10 * .e14
#   .e345 <- .e10 * .e256 - .e281 * .e6 * .e7/.e1
#   .e348 <- .e62 * .e11/.e16 + y * .e272 * .e60/.e1
#   .e349 <- .e5^(.e4 + .e18)
#   .e350 <- 2 * .e15
#   .e351 <- 2 * .e17
#   .e352 <- .e314 + .e27
#   .e353 <- xi * .e105
#   .e354 <- .e283 * .e62
#   .e355 <- xi * .e30
#   .e356 <- .e317 * .e57
#   .e357 <- .e319 * .e57
#   .e358 <- .e235 * .e57
#   .e359 <- .e235 * .e38
#   .e360 <- .e321 * .e57
#   .e362 <- .e237 * .e57 + .e294 * .e33/.e61
#   .e363 <- .e325 * .e57
#   .e364 <- .e241 * .e57
#   .e365 <- .e170 * .e93
#   .e366 <- .e284 - .e259
#   .e367 <- .e323 - .e328
#   .e368 <- ((.e19 + .e196/.e39) * .e6 - .e19/xi) * .e62
#   .e370 <- .e244 * .e57 + .e290 * .e33/.e61
#   .e371 <- .e188 * .e93
#   .e372 <- .e285 - .e245 * .e24
#   .e373 <- .e188 * .e38
#   .e374 <- .e260 * .e19
#   .e375 <- .e286 - .e307/.e25
#   .e376 <- .e247 * .e46
#   .e377 <- .e183 * .e24
#   .e378 <- .e287 - .e300/.e25
#   .e379 <- (((.e95 - .e249)/xi + .e36) * .e11 + .e254)/xi
#   .e380 <- .e191 * .e93
#   .e381 <- .e288 - .e242 * .e24
#   .e382 <- .e191 * .e38
#   .e383 <- .e290/.e25
#   .e384 <- .e242 * .e98
#   .e385 <- .e242 * .e81
#   .e386 <- .e292/.e10
#   .e387 <- .e259 * .e81
#   .e388 <- .e294/.e25
#   .e389 <- .e245 * .e98
#   .e390 <- .e245 * .e81
#   .e391 <- .e146 * .e38
#   .e393 <- .e302 * .e12 + .e59/.e10
#   .e395 <- .e304 * .e12 + .e59/.e39
#   .e396 <- .e304 * .e33
#   .e397 <- .e306 * .e48
#   .e398 <- (.e339 + y * ((.e10 + .e343 * .e11)/xi - (.e10 + .e315) * .e14)) * .e6
#   .e399 <- .e165^2
#   .e400 <- .e309 * .e57
#   .e401 <- .e345 * .e57
#   .e406 <- .e109 * .e46 * .e13
#   .e410 <- (.e17 - 2 * (.e78^2 * .e17/.e56)) * .e1 + .e2 * .e14 * (.e315 - .e10)
#   .e415 <- .e57 * .e227
#   .e416 <- .e57 * .e256
#   .e417 <- .e19 * .e119
#   .e418 <- .e101 * .e163
#   .e419 <- .e15 * .e7
#   .e424 <- .e97 * .e12
#   .e425 <- .e97 * .e24
#   .e426 <- .e337/.e312
#   .e427 <- .e41 * .e38
#   .e428 <- .e58 + .e350
#   .e430 <- .e41/.e264 + .e106 * .e272 * .e60
#   .e431 <- .e296 + .e300 * .e33/.e61
#   .e434 <- .e299 * .e33/.e61 - .e295
#   .e437 <- .e141/.e61 - 2 * (.e15 * .e30 * .e119 * .e267/.e92)
#   .e438 <- .e305 + .e307 * .e33/.e61
#   .e439 <- .e46 * .e1
#   .e440 <- .e262 * .e33
#   .e442 <- (1 - 2 * (.e104 * .e10 * .e1/.e64)) * .e1 * .e38
#   .e452 <- .e5^(.e112 - 1)
#   .e454 <- .e39 + .e280 - .e279
#   .e455 <- .e5^.e100
#   .e458 <- .e227 * .e38
#   .e461 <- .e38 * .e256
#   .e462 <- .e33^2
#   .e463 <- .e115^2
#   .e464 <- .e167^2
#   .e465 <- .e135^2
#   .e466 <- .e168^2
#   .e467 <- .e144^2
#   .e468 <- 1/.e5
#   .e469 <- 2 * (.e171 * .e30 * .e1/.e92)
#   .e471 <- y * (.e109 * .e14 - .e10/xi)
#   .e473 <- y * (.e6 * (.e352 - .e29) - .e250)/.e1
#   .e476 <- y * (.e14 * .e38 - .e10/.e16)/.e1 - ((.e11 * (.e74 - (.e73/xi + .e351)/xi) + .e71/.e1)/.e16 + 2 * (.e17 * .e267 * .e77^2/.e56))
#   .e480 <- y * .e310 * .e63/.e1 - 2 * (.e69 * .e11/.e16)
#   .e483 <- .e2 * .e310 * .e63/.e1
#   .e484 <- y/(.e5 * .e1)
#   .e485 <- xi * .e237
#   .e486 <- xi * .e241
#   .e489 <- xi * .e170 * .e6 * .e7
#   .e490 <- xi * .e244
#   .e494 <- xi * .e188 * .e6 * .e7/.e1
#   .e497 <- .e173 * .e15 * .e6 - .e299/.e25
#   .e501 <- xi * .e191 * .e6 * .e7/.e1
#   .e503 <- xi * .e109 * .e15
#   .e506 <- xi * .e348 * .e19 * .e6
#   .e507 <- xi * .e77
#   .e508 <- .e16 * .e6
#   .e510 <- .e134 * .e15 * .e6
#   # final matrix
#   out <- matrix(0, length(y), 25)
#   # third derivatives
#   out[, 1] <- y * ((.e10 * (xi * ((.e410/.e56 - 
#         (.e62 + .e2 * .e272 * .e60/.e1)/.e1) * .e19 - y * .e137 * 
#         .e15 * .e12/.e25) * .e6 - (((.e50 + .e42 * (2 * .e483 - 
#         2 * .e69)) * .e15 - .e147 * (.e469 + .e202/.e25))/.e92 + 
#         (y * (.e497/.e10 + .e354 * .e53/.e1 - ((.e46/.e142 + 
#             y * .e257 * .e44/(.e5^(2 + .e18 + .e31) * .e1))/.e1 + 
#             (.e165 - .e2 * .e132 * .e349) * .e46/.e399) * .e13) - 
#             .e178)/.e1) * .e12) - .e107 * (y * (2 * (.e7 * .e256) - 
#         .e84 * .e22/.e1) - .e84 * .e7)/.e1)/.e19 + .e12 * (y * 
#         (.e401/.e10 + .e416 - .e383) - .e242)/.e1) #d111
#   out[, 2] <- y * 
#         ((((((((.e10 - 2 * (.e104^2 * .e10/.e64)) * .e1 + .e42 * 
#             (.e27 - .e7))/.e64 + (.e49 + xi * (y * (xi * ((2 * 
#             (xi * .e5^(.e100 - 2)/.e67) - .e108 * .e20) * .e11 - 
#             .e108)/.e67 - .e6 * .e62)/.e1 - .e116/.e67))/.e1) * 
#             .e15 - .e13 * (y * (.e375/.e10 + .e286/.e10)/.e1 - 
#             .e104 * .e46 * .e48/.e64))/.e43 + (y * (.e378/.e10 + 
#             .e287/.e10 + xi * (.e15 * .e62 * .e55/.e1 + xi * 
#             (.e152 * .e7 + (2 * (.e510 * .e123/.e82) - (.e15 * 
#                 .e22 + .e59/.e5)) * .e48/.e1)/.e82) * .e6) + 
#             .e508 * (y * .e174 * .e7 - .e419 * .e48)/.e82 - .e180)/.e1) * 
#             .e12 - (.e410 * .e6/.e56 + (.e149 + .e16 * (y * (xi * 
#             ((2 * (xi^5 * .e6 * .e123/.e79) - .e22) * .e11 - 
#                 .e96) * .e6/.e79 - .e343/.e465)/.e1 - .e107 * 
#             .e7 * .e11/.e79))/.e1) * .e19) * .e10 - .e107 * (y * 
#             (2 * (.e235 * .e7) - .e90 * .e22/.e1) - .e90 * .e7)/.e1)/.e19 + 
#             .e12 * (y * (.e356/.e10 + .e358 - .e388) - .e245)/.e1) #d112
#   out[, 3] <- -(y * (y * ((.e400/.e10 + .e415 - .e386) * .e12 - 
#             ((((.e12 * (.e354 * .e24 - .e302/.e10) + .e354 - 
#                 .e59/.e30)/.e10 + xi * (.e146 * .e62 + .e430 * 
#                 .e24 + .e15/.e264) * .e6) * .e8 + (2 * (.e283/.e264) - 
#                 .e59/.e142) * .e12 + xi * .e430 * .e6) * .e10 + 
#                 .e107 * (2 * (.e7 * .e227) - .e87 * .e22))/.e19)/.e1 - 
#             (.e259 * .e12 + .e309/.e19))/.e1) #d113
#   out[, 4] <- y * 
#         ((((((((.e375 * .e33 + .e376 * .e48/.e1)/.e61 + (2 * 
#             (xi * .e46 * .e251 * .e48 * .e33/.e313) - .e305/.e10)/.e1) * 
#             .e13 + ((.e442 + y * .e260)/.e64 - (((.e216 - (2 * 
#             .e22 + .e233)) * .e11 + .e353 * (2 * (xi * .e251 * 
#             .e11/.e67) - .e468))/(.e1 * .e67) + xi * .e104/.e463)) * 
#             .e15)/.e10 + (xi * .e172 * .e6 * .e62 - 2 * (.e198/.e265))/.e1 - 
#             (.e378 * .e33 + .e326 * .e55/.e1)/.e39)/xi + .e296/.e25 - 
#             ((.e152 + 2 * (.e510 * .e265 * .e48/(.e1 * .e82))) * 
#                 .e109 + .e48 * (xi * (.e246 - .e6 * .e352) * 
#                 .e15 - .e406/.e10)/.e1)/.e82) * .e12 - (.e398/.e56 - 
#             .e16 * (((.e246 - (.e7 + .e314 + .e27) * .e6) * .e11 + 
#                 .e133 * (2 * (xi^4 * .e6 * .e265 * .e11/.e79) - 
#                   .e468))/(.e1 * .e79) + 2 * (.e78/.e467))) * 
#             .e19) * .e10 + .e359 + (.e260 * .e90 - .e485 * .e6 * 
#             .e7)/.e1)/.e19 + (.e357/.e25 - (((.e358 - .e388) * 
#             .e98 - .e282 * .e57 * .e83/.e276) * .e33 + .e247 * 
#             .e90 * .e57 * .e98/.e1)/xi) * .e12) #d122
#   out[, 5] <- -(y * 
#         (((((.e247 * .e139 + ((2 * (.e297/.e265) - .e302/.e39) * 
#             .e12 + 2 * (.e15/.e265) - .e59/.e123) * .e33)/.e43 + 
#             .e107 * (xi * .e139 * .e7 * .e33/.e82 - 2 * (.e192/.e264)) - 
#             (.e374 * .e24 - .e391)/.e30) * .e8 + ((.e326 - .e141/.e10)/.e167 + 
#             .e182 * .e111 * .e455 * .e33/.e464) * .e12 - ((.e374 - 
#             .e427/.e10)/.e30 + 2 * (.e114 * .e38/.e264))) * .e10 + 
#             .e260 * .e87 + .e458 - .e486 * .e6 * .e7)/.e19 + 
#             (.e360/.e10 - (((.e415 - .e386) * .e98 - .e273 * 
#                 .e57 * .e83/.e7) * .e33 + .e247 * .e87 * .e57 * 
#                 .e98)/xi) * .e12)/.e1) #d123
#   out[, 6] <- -(y * 
#             ((.e284 * .e12 + (.e10 * (xi * (.e156 * .e62 + .e19 * 
#                 .e62) * .e6 - ((.e146 * .e8 + .e41/.e10)/.e10 + 
#                 .e197)) - .e489)/.e19 + .e273 * .e93 * .e6 * 
#                 .e7 - (.e93 * .e227 - .e387/.e10) * .e10) * .e24 - 
#                 .e367/(.e9 * .e10)) * .e8/.e1) #d133
#   out[, 7] <- (((((((.e438 * .e33 + .e46 * .e48 * 
#         .e255)/.e61 + (.e305/.e61 - .e454 * .e46 * .e48/.e313) * 
#         .e33) * .e13 + .e15 * (y * ((.e270/xi - 2 * (.e10 * .e267 * 
#         .e38^2/.e64))/.e64 - .e109/.e463) * .e1 - (.e105 * (.e484 - 
#         2 * (.e353 * .e7 * .e11/.e67)) + .e11 * (.e231 - (.e91 - 
#         .e224) * .e11/.e16))/.e67))/.e10 + ((.e426 + .e130) * 
#         .e55 * .e33 - .e185)/xi - ((.e431 * .e33 + .e102 * .e255)/.e39 + 
#         .e172 * .e38/.e30 + .e296 * .e33/.e39))/xi - ((.e406 * 
#         .e33/.e61 + .e15 * (.e473 - .e379)) * .e48 + (.e164 - 
#         2 * (.e503 * .e10 * .e48/.e82)) * .e109)/.e82) * .e12 - 
#         .e19 * (y * ((.e6 * .e476 - .e77/.e16)/.e56 - 2 * (xi * 
#             (.e351 + .e507)/.e467)) * .e1 - ((.e95 + xi * (2 * 
#             .e38 + .e473 - .e379)) * .e11 + .e151 * (.e484 - 
#             2 * (.e134 * .e10 * .e133 * .e11/.e79)))/.e79)) * 
#         .e10 + .e90 * .e270/xi + 2 * (.e237 * .e38))/.e19 - ((.e362 * 
#         .e98 + .e357/.e39 + .e90 * .e342 * .e57) * .e33 + .e336 * 
#         .e90) * .e12/xi #d222
#   out[, 8] <- -((((((((.e396/.e39 - .e337 * 
#         .e24/.e312) * .e12 + .e141/.e312 - .e426) * .e33/xi + 
#         .e139 * .e255)/.e10 - (.e139 * .e38 * .e33 + .e101 * 
#         .e270)/.e30)/xi + .e417 * .e24 * .e38/.e311 - .e139 * 
#         .e109 * .e33/.e82) * .e8 + ((.e15 * .e255 + .e59 * .e462/.e61)/.e167 - 
#         (.e123 + xi * (y * .e111 * .e455/.e1 - 3 * (.e123 * .e11/.e16))) * 
#             .e15 * .e33/.e464) * .e12 + .e417 * .e38/.e311 - 
#         (.e19 * .e270 + .e427 * .e33/.e39)/.e355) * .e10 + .e87 * 
#         .e270/xi + 2 * (.e241 * .e38))/.e19 - (((.e364 + .e292 * 
#         .e33/.e61) * .e98 + .e360/.e39 + .e87 * .e342 * .e57) * 
#         .e33 + .e336 * .e87) * .e12/xi) #d223
#   out[, 9] <- -(((((((.e139 * 
#         .e8 + .e41/.e39)/.e10 + .e41/.e123) * .e33/xi - (.e156/.e30 + 
#         .e19/.e30) * .e38) * .e10 + .e170 * .e38)/.e19 - ((.e241 * 
#         .e93 + .e387 * .e33/.e61) * .e10 + .e284 * .e98 * .e12 * 
#         .e33/xi + .e275 * .e38)) * .e24 + .e367 * .e33/(xi * 
#         .e9 * .e39)) * .e8) #d233
#   out[, 10] <- -(((((.e66 + 2 * .e19)/.e10 + 
#         .e68) * .e10/.e19 - ((.e365 - .e275) * .e10 + .e365 * 
#         .e10)) * .e8 * .e24 + .e323 - .e328) * .e8 * .e24) #d333
#   # fourth derivatives
#   ## ignored for now
#   out
# }
    
.egpd1.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd1d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd1.d12 <- function(pars, likdata) {
egpd1d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd1.d34 <- function(pars, likdata) {
egpd1d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.iG1 <- function(v, kappa) v^(1/kappa)

.egpd1fns <- list(d0=.egpd1.d0, d120=.egpd1.d12, d340=.egpd1.d34, m=1, iG=.iG1)

## model 2 ##

.G2 <- function(v, kappa1, kappa2, p) {
p * v^kappa1 + (1 - p) * v^kappa2
}

.iG2i <- function(v, kappa1, kappa2, p) {
vv <- range(c(v^1/kappa1, v^1/kappa2))
d <- diff(vv)
lo <- vv[1]
while(.G2(lo, kappa1, kappa2, p) - v > 0) lo <- max(0, lo - d)
hi <- vv[2]
while(.G2(hi, kappa1, kappa2, p) - v > 0) hi <- min(1, hi + d)
uniroot(function(x) .G2(x, kappa1, kappa2, p) - v, c(lo, hi))$root
}

.iG2 <- function(v, kappa1, kappa2, p) {
n <- length(v)
vapply(seq_len(n), function(i) .iG2i(v[i], kappa1[i], kappa2[i], p[i]), double(1))
}

.egpd2.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd2d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd2.d12 <- function(pars, likdata) {
egpd2d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd2.d34 <- function(pars, likdata) {
egpd2d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd2fns <- list(d0=.egpd2.d0, d120=.egpd2.d12, d340=NULL, m=2, iG=.iG2)

## model 3 ##

.egpd3.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd3d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd3.d12 <- function(pars, likdata) {
egpd3d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd3.d34 <- function(pars, likdata) {
egpd3d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.iG3 <- function(v, delta) 1 - qbeta(1 -v, 1/delta, 2)^(1/delta)

.egpd3fns <- list(d0=.egpd3.d0, d120=.egpd3.d12, d340=.egpd3.d34, m=3, iG=.iG3)

## model 4 ##

.egpd4.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd4d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd4.d12 <- function(pars, likdata) {
egpd4d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd4.d34 <- function(pars, likdata) {
egpd4d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.iG4 <- function(v, delta, kappa) 1 - qbeta(1 - v^(2/kappa), 1/delta, 2)^(1/delta)

.egpd4fns <- list(d0=.egpd4.d0, d120=.egpd4.d12, d340=NULL, m=4, iG=.iG4)

