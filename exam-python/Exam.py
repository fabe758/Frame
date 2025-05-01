from setup_o import *
# from setup_s import *

ml41 = RT.dMlens({RT.dLens(0.9, RT.dV2d(0.0, 0.0)),
                  RT.dLens(0.04, RT.dV2d(1.13, 0.11)),
                  RT.dLens(0.04, RT.dV2d(0.98, -0.21)),
                  RT.dLens(0.02, RT.dV2d(1.22, -0.22))});

src01 = RT.dSource(RT.dCircle(RT.dV2d(0.1, 0.0), 0.01))
fr41 = RT.dFrTri(ml41, src01)
t0 = RT.dTime(0.0)

ln1 = RT.dLinear(-0.3, t0, 4.0 * RT.atan(1.0) * 2.0, 20.0)

sm01 = RT.dSrcMotion(src01, ln1)
sm1 = RT.dSrcMotion(src01, RT.dLinear())

# mlm41 = RT.MlMotion(ml41)

# lc41_01 = RT.dLcurve({sm01}, RT.MlMotion(ml41), t.range(200, -20.0, 20.0))
# lc41_01 = RT.dLcurve({sm01}, mlm41, t.range(200, -20.0, 20.0))
lc41_01 = RT.dLcurve({sm01}, ml41, t0.range(200, -20.0, 20.0))

lc41_01.process_while(1000)

rtlc = RT.RTLcurve(lc41_01)

mag = rtlc.mag()

# mag.Draw()
#
# aa = input()
