
import system/ansi_c

type
  Mpfr {.importc: "__mpfr_struct", header: "mpfr.h".} = object
  mpfr_prec_t = cint
  mpfr_ptr = ptr Mpfr
  mpfr_srcptr = ptr Mpfr

  mpfr_rnd_t = enum
    MPFR_RNDNA = -1 # round to nearest, with ties away from zero (mpfr_round)
    MPFR_RNDN = 0,  # round to nearest, with ties to even 
    MPFR_RNDZ,      # round toward zero 
    MPFR_RNDU,      # round toward +Inf
    MPFR_RNDD,      # round toward -Inf 
    MPFR_RNDA,      # round away from zero
    MPFR_RNDF,      # faithful rounding


proc mpfr_init2(m: mpfr_ptr, prec: mpfr_prec_t) {.importc, header: "mpfr.h".}
proc mpfr_set_ui(m: mpfr_ptr, val: culong, rnd: mpfr_rnd_t): cint {.importc, header: "mpfr.h".}
proc mpfr_add(m: mpfr_ptr, v1, v2: mpfr_ptr, rnd: mpfr_rnd_t): cint {.importc, header: "mpfr.h".}
proc mpfr_set_si(m: mpfr_ptr, val: clong, rnd: mpfr_rnd_t): cint {.importc, header: "mpfr.h".}
proc mpfr_pow(m: mpfr_ptr, s: mpfr_srcptr, v: mpfr_ptr, rnd: mpfr_rnd_t): cint {.importc, header: "mpfr.h".}
proc mpfr_pow_ui(m: mpfr_ptr, s: mpfr_srcptr, v: culong, rnd: mpfr_rnd_t): cint {.importc, header: "mpfr.h".}
proc mpfr_set_d(m: mpfr_ptr, v: cdouble, rnd: mpfr_rnd_t): cint {.importc, header: "mpfr.h".}
proc mpfr_sin(m: mpfr_ptr, s: mpfr_srcptr, rnd: mpfr_rnd_t) {.importc, header: "mpfr.h".}
proc mpfr_printf(fmt: cstring) {.importc, header: "mpfr.h", varargs.}
proc mpfr_asprintf(p: ptr cstring, fmt: cstring) {.importc, header: "mpfr.h", varargs.}
proc mpfr_clear(m: mpfr_ptr) {.importc, header: "mpfr.h".};



var ctxPrec = 53.mpfr_prec_t
var ctxRnd = MPFR_RNDN

proc initMpfr(): Mpfr =
  mpfr_init2(result.addr, ctxPrec)

converter toMpfr(v: float): Mpfr =
  result = initMpfr()
  discard mpfr_set_d(result.addr, v, ctxRnd)

converter toMpfr(v: int): Mpfr =
  result = initMpfr()
  discard mpfr_set_si(result.addr, v, ctxRnd)

proc `+`(m1, m2: Mpfr): Mpfr =
  result = initMpfr()
  discard mpfr_add(result.addr, m1.unsafeAddr, m2.unsafeAddr, ctxRnd)

proc `^`(m: Mpfr, v: int): Mpfr =
  result = initMpfr()
  discard mpfr_pow_ui(result.addr, m.unsafeAddr, v.culong, ctxRnd)

proc `^`(m: Mpfr, v: float): Mpfr =
  result = initMpfr()
  var v2 = Mpfr(v)
  discard mpfr_pow(result.addr, m.unsafeAddr, v2.addr, ctxRnd)

proc sin(m: Mpfr): Mpfr =
  result = initMpfr()
  mpfr_sin(result.unsafeAddr, m.unsafeAddr, ctxRnd)

proc `$`(m: var Mpfr): string =
  var buf: cstring
  mpfr_asprintf(buf.addr, "%.50Rg", m.addr)
  result = $buf
  c_free(buf)



let y = 5.Mpfr
var x = sin((y + y) ^ 22.0)
echo $x


