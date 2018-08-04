#!/bin/env python

"""
Extract the spectrum and OPE coefficients from an sdp file and the
output of sdpb.  Run './spectrum.py --help' for usage information.
This script will only give sensible output when the output file
contains a primal-dual optimal solution.

Requirements:

  - The arbitrary precision polynomial root-finder MPSolve:

    http://numpi.dm.unipi.it/mpsolve-2.2/

    The path to the MPSolve executable 'unisolve' can be passed as an
    argument to this script.
    
  - mpmath for arbitrary precision arithmetic:

    http://mpmath.org/

This script was originally written for use in:

  - Komargodski, Zohar and Simmons-Duffin, David, 'The Random Bond
    Ising Model in 2.01 and 3 Dimensions,' arXiv:1603.04444

An explanation of how it works appears in:

  - Simmons-Duffin, David, 'The Lightcone Bootstrap and the Spectrum
    of the 3d Ising CFT,' arXiv:1612.08471

Author(s): David Simmons-Duffin

Copyright 2016, David Simmons-Duffin.  Distributed under the MIT
license.  (See accompanying file LICENSE or copy at
http://opensource.org/licenses/MIT)
"""

from mpmath import *
import argparse
import itertools
import lxml.etree
import operator
import os
import subprocess
import sys
import tempfile

ZERO_ROOT_THRESHOLD = 1e-30
"""Float: threshold for deciding whether a polynomial has a root at
zero. Specifically, we test whether p(0)/p'(0) < ZERO_ROOT_THRESHOLD.
"""

UNISOLVE_EXECUTABLE = os.environ["HOME"]+"/bin/unisolve"
"""Path to the unisolve executable
"""

################# parsing SDPB output files #################

def parseOutFile(outFile):
    """Read outFile, assumed to be in the SDPB output format described
    in the SDPB manual"""
    def readArray(str):
        return map(mpf, str.rstrip("}").lstrip("{").split(", "))
    parseFields = {
        "terminateReason": lambda x: x.strip("\""),
        "primalObjective": mpf,
        "dualObjective": mpf,
        "dualityGap": mpf,
        "primalError": mpf,
        "dualError": mpf,
        "runtime": float,
        "y": readArray,
        "x": readArray
        }
    output = {}
    for l in open(outFile):
        [key, val] = map(lambda s: s.strip(" ;\n"), l.split(" = "))
        output[key] = parseFields[key](val)
    return output

################# parsing SDP xml files #################

# This section implements a streaming parser for
# PolynomialVectorMatrices. The main point is to avoid reading the
# entire xml file into memory at once. Instead, a single
# PolynomialVectorMatrix can be extracted and processed at a time.

class Polynomial(object):
    """A basic Polynomial class. Should perhaps be replaced with an
    implementation from a library.
    """
    def __init__(self, coeffs):
        self.coeffs = coeffs
    def __add__(self, other):
        if isinstance(other, Polynomial):
            return Polynomial([
                a + b for a, b in itertools.izip_longest(self.coeffs, other.coeffs, fillvalue=0)
                ])
        else:
            return self + Polynomial([other])
    def __neg__(self):
        return self*(-1)
    def __sub__(self, other):
        return self + -other
    def __rsub__(self, other):
        return -self + other
    def __mul__(self, other):
        if isinstance(other, Polynomial):
            coeffs = [0]*(len(self.coeffs)+len(other.coeffs) - 1)
            for i1, c1 in enumerate(self.coeffs):
                for i2, c2 in enumerate(other.coeffs):
                    coeffs[i1+i2] += c1*c2
            return Polynomial(coeffs)
        else:
            return Polynomial([c * other for c in self.coeffs])
    def __str__(self):
        return " + ".join(str(c) + " x^" + str(i) for i, c in enumerate(self.coeffs))
    def __repr__(self):
        return "Polynomial(" + repr(self.coeffs) + ")"
    def val(self, x):
        return reduce(lambda acc, c: acc * x + c, reversed(self.coeffs), 0)
    def degree(self):
        return len(self.coeffs) - 1
    def derivative(self):
        coeffs = [c*i for i, c in enumerate(self.coeffs)]
        return Polynomial(coeffs[1:])

class PolynomialVectorMatrix(object):
    def __init__(self, rows, cols, elements, samplePoints, sampleScalings, bilinearBasis):
        self.rows           = rows
        self.cols           = cols
        self.elements       = elements
        self.samplePoints   = samplePoints
        self.sampleScalings = sampleScalings
        self.bilinearBasis  = bilinearBasis

class ParseException(Exception):
    """Used if a function encounters a tag it doesn't understand. By
    raising ParseException, it can pass the tag up the call stack to a
    parent routine that knows how to handle the tag.
    """
    def __init__(self, elem):
        self.elem = elem

def checkAndClear(stream, tag, f):
    """When the parser encounters a close tag, check that the tag name
    is the expected one. If so, compute f(elem) and clear elem from
    the parse tree (freeing up memory). If not, throw ParseException
    with the unexpected tag to pass control to another function
    """
    (event, elem) = stream.next()
    if event == "end" and elem.tag == tag:
        result = f(elem)
        elem.clear()
        return result
    else:
        raise ParseException(elem)

def parseElem(stream, tag):
    return checkAndClear(stream, tag, lambda e: e.text)

def parseArray(stream, tag, parseElt):
    a = []
    while True:
        try:
            x = parseElt(stream)
            a.append(x)
        except ParseException as e:
            if e.elem.tag == tag:
                e.elem.clear()
                return a
            else:
                raise e

def parseReal(stream, tag):
    return mpf(parseElem(stream, tag))

def parseInt(stream, tag):
    return int(parseElem(stream, tag))

def parsePolynomial(stream):
    return Polynomial(parseArray(stream, "polynomial", lambda s: parseReal(s, "coeff")))

def parsePolynomialVector(stream):
    return parseArray(stream, "polynomialVector", parsePolynomial)

def parsePolynomialVectorMatrix(stream):
    m = PolynomialVectorMatrix(parseInt(stream, "rows"),
                               parseInt(stream, "cols"),
                               parseArray(stream, "elements", parsePolynomialVector),
                               parseArray(stream, "samplePoints", lambda s: parseReal(s, "elt")),
                               parseArray(stream, "sampleScalings", lambda s: parseReal(s, "elt")),
                               parseArray(stream, "bilinearBasis", parsePolynomial))
    return checkAndClear(stream, "polynomialVectorMatrix", lambda e: m)

# We won't use this because it would read the entire SDP into memory,
# but this is how one can parse the whole SDP at once
class SDP(object):
    def __init__(self, stream):
        self.objective                = parseArray(stream, "objective", lambda s: parseReal(s, "elt"))
        self.polynomialVectorMatrices = parseArray(stream, "polynomialVectorMatrices", parsePolynomialVectorMatrix)

################# computation helpers #################

def dotPolVector(polVector, alpha):
    return sum((p * a for p, a in zip(polVector, alpha)), Polynomial([]))

def partition(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def prod(iterable):
    return reduce(operator.mul, iterable, 1)

def indexMax(values):
        return max(xrange(len(values)),key=values.__getitem__)

def lagrangeInterpolationMatrix(xs, ys):
    return matrix([
        [ prod(((y - xk)/(xi - xk) if k != i else 1)
               for k, xk in enumerate(xs))
          for y in ys
          ]
        for i, xi in enumerate(xs)
        ])

def det(m):
    n = len(m)
    if n == 1:
        return m[0][0]
    else:
        result = 0
        for j in range(n):
            result += ((-1)**j) * m[j][0] * det([row[1:] for row in m[:j] + m[j+1:]])
        return result

def symmetricMatrix(xs, rows):
    """Construct an upper-triangular matrix from a list of entries, and symmetrize it
    """
    m = mp.matrix(rows, rows)
    xIter = iter(xs)
    for s in range(rows):
        for r in range(s+1):
            m[r,s] = xIter.next() / (2 if r != s else 1)
            m[s,r] = m[r,s]
    return m

def upperTriangularPart(m):
    """The inverse of 'symmetricMatrix'
    """
    for s in range(m.rows):
        for r in range(s+1):
            yield m[r,s] * (2 if r != s else 1)

def sqrtRankOneComponent(m):
    """Find v such that the rank-1 approximation m = v v^T is best. If
    the largest eigenvalue of m is negative, return None.
    """
    e, q = mp.eigsy(m)
    i = indexMax(e)
    if e[i] > 0:
        return sqrt(e[i])*q[:,i]
    else:
        return None

def outerProduct(xs,ys):
    m = mp.matrix(len(xs), len(ys))
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            m[i,j] = x*y
    return m

def leastSquaresFit(m):
    """Return a function taking a vector v to the least-squares
    solution of m x = v
    """
    U, S, V = mp.svd_r(m)
    SInv = mp.diag(map(lambda x: 1/x, S))
    return (lambda v: V.T*(SInv*(U.T*v)))

################# polynomial root finding  #################

def writePolFile(p, f):
    """Write a polynomial to a file in the format required for
    unisolve
    """
    f.write("drf\n")
    f.write(str(mp.dps) + "\n")
    f.write(str(len(p.coeffs) - 1) + "\n")
    for c in p.coeffs:
        f.write(str(c) + "\n")
    f.flush()

def unisolveRoots(p):
    """Call UNISOLVE_EXECUTABLE to determine the roots of p
    """
    with tempfile.NamedTemporaryFile(suffix=".pol") as f:
        writePolFile(p, f)
        out = subprocess.check_output([
            UNISOLVE_EXECUTABLE,
            "-H1",
            "-o" + str(mp.dps),
            "-Oc",
            "-Ga",
            f.name])
        return [mpc(*n.strip("()").split(", "))
                for n in out.rstrip('\n').split('\n')]

def posRealRoots(p):
    return [r.real for r in unisolveRoots(p)
            if r.imag == 0 and r.real > 0]

################# extracting operators  #################

def determinantPositiveRoots(m, zeroRootThreshold):
    """Given a polynomial matrix m, positive semidefinite for x >= 0,
    estimate the roots of m. Since p = det(m) never actually hits
    zero, we estimate the roots by computing the local minima. We then
    check whether p(x)/p''(x) < zeroRootThreshold. To check for a
    possible additional root at zero, we check whether p(0)/p'(0) <
    zeroRootThreshold.
    """
    d = det(m)
    if d.degree() == 0:
        return [mpf(0)]
    else:
        dPrime = d.derivative()
        dPrime2 = dPrime.derivative()
        rs = [r for r in posRealRoots(dPrime)
              if (dPrime2.val(r) >= 0 and
                  d.val(r)/dPrime2.val(r) < zeroRootThreshold)
              ]
        if abs(d.val(0)/dPrime.val(0)) < zeroRootThreshold:
            rs = [mpf(0)] + rs
        return rs

def operatorsFromPolVectorMatrix(m, alpha, xStream):
    mDotAlpha = partition([dotPolVector(v, alpha) for v in m.elements], m.rows)
    xVectors = [
        [s*x for s,x in zip(m.sampleScalings, xStream)]
        for i in range((m.rows*(m.rows+1))/2)
        ]
    operators = []
    error = mp.matrix(xVectors)
    roots = determinantPositiveRoots(mDotAlpha, ZERO_ROOT_THRESHOLD)
    if roots:
        rootsInterpMat = lagrangeInterpolationMatrix(m.samplePoints, roots)
        rootsFit = leastSquaresFit(rootsInterpMat)
        lambdaEntries = [rootsFit(mp.matrix(xs)) for xs in xVectors]
        for i, r in enumerate(roots):
            lambdaMat = symmetricMatrix(map(lambda v: v[i], lambdaEntries), m.rows)
            lambdaVec = sqrtRankOneComponent(lambdaMat)
            if lambdaVec is not None:
                operators.append((r, lambdaVec))
                error -= outerProduct(
                    list(upperTriangularPart(outerProduct(lambdaVec,lambdaVec))),
                    rootsInterpMat[:,i])
    return (operators, mp.norm(error))

def readSpectrum(sdpStream, output):
    """Read a PolynomialVectorMatrix one at a time from sdpStream and
    extract the associated operators, using the functional and primal
    solution from output.
    """
    alpha = [1] + output["y"]
    xStream = iter(output["x"])
    # we don't actually use the objective, but we must parse it from
    # sdpStream first before parsing the polynomialVectorMatrices
    _ = parseArray(sdpStream, "objective", lambda s: parseReal(s, "elt"))
    while True:
        try:
            m = parsePolynomialVectorMatrix(sdpStream)
            yield operatorsFromPolVectorMatrix(m, alpha, xStream)
        except ParseException as e:
            if e.elem.tag != "polynomialVectorMatrices":
                raise e
            else:
                break

def toMathReal(x):
    """mpfloat to Mathematica string."""
    return str(x).replace('e', '*^')

def toMathArray(xs, toMathElt):
    """Array to mathematica string."""
    return '{' + ", ".join(toMathElt(x) for x in xs) + '}'

def toMathOp(op):
    """Operator to mathematica string."""
    delta, lamdaVec = op
    return '{' + toMathReal(delta) + ', ' + toMathArray(lamdaVec, toMathReal) + '}'

def main():
    global ZERO_ROOT_THRESHOLD
    global UNISOLVE_EXECUTABLE
    argParser = argparse.ArgumentParser(
        description=('Extract spectrum and OPE coefficients from an sdp file '
                     'and an sdpb output file. Uses the arbitrary precision '
                     'polynomial solver \'unisolve\'.'))
    argParser.add_argument('-s', '--sdpFile', required=True,
                           help='sdp file (.xml format)')
    argParser.add_argument('-o', '--outFile',
                           help='sdpb output file (.out format)')
    argParser.add_argument('-m', '--spectrumFile',
                           help='spectrum output file (.m format)')
    argParser.add_argument('-p', '--precision', type=int, default=200,
                           help='working precision in decimal digits (default 200)')
    argParser.add_argument('-z', '--zeroRootThresh', type=float,
                           default=ZERO_ROOT_THRESHOLD,
                           help=('threshold for p(0)/p\'(0) to determine if p '
                                 'has a root at zero '
                                 '(default '+str(ZERO_ROOT_THRESHOLD)+')'))
    argParser.add_argument('-u', '--unisolveExec', default=UNISOLVE_EXECUTABLE,
                           help=('path to the unisolve executable '
                                 '(default '+UNISOLVE_EXECUTABLE+')'))
    args = argParser.parse_args()
    if args.outFile is None:
        args.outFile = args.sdpFile.replace(".dat-s", ".out").replace(".xml", ".out")
    if args.spectrumFile is None:
        args.spectrumFile = args.outFile.replace(".out", ".spectrum.m")

    mp.dps = args.precision
    ZERO_ROOT_THRESHOLD = args.zeroRootThresh
    UNISOLVE_EXECUTABLE = args.unisolveExec

    print "reading sdp file     :", args.sdpFile
    sdpStream = lxml.etree.iterparse(args.sdpFile)

    print "reading out file     :", args.outFile
    output = parseOutFile(args.outFile)
    if output["terminateReason"] != 'found primal-dual optimal solution':
        raise Exception("terminateReason should be 'found primal-dual "
                        "optimal solution', not '" + output["terminateReason"] + "'")

    print "writing spectrum file:", args.spectrumFile
    sys.stdout.flush()
    spec = open(args.spectrumFile, 'w')
    spec.write("{")
    for j, (ops, err) in enumerate(readSpectrum(sdpStream, output)):
        print "---- j =", j
        print "deltas =", [float(x) for x, _ in ops]
        print "error  =", float(err), "(" + str(len(ops)) + " operators)"
        sys.stdout.flush()
        spec.write(str(j) + " -> " + toMathArray(ops, toMathOp))
        spec.write(",\n")
        spec.flush()
    spec.write('"objective" -> ' + toMathReal(output['primalObjective']))
    spec.write("}")
    spec.close()

if __name__ == '__main__':
    main()
