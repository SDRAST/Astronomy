# -*- coding: utf-8 -*-
"""
from http://stackoverflow.com/questions/2417794/how-to-make-the-angles-in-a-matplotlib-polar-plot-go-clockwise-with-0-at-the-top
"""

import numpy as N
import matplotlib.pyplot as P

import matplotlib.projections as MPLproj
import matplotlib.transforms as MPLtrans

class NorthPolarAxes(MPLproj.PolarAxes):
    '''
    A variant of PolarAxes where theta starts pointing north and goes
    clockwise.
    '''
    name = 'northpolar'

    class NorthPolarTransform(MPLproj.PolarAxes.PolarTransform):
        def transform(self, tr):
            xy   = N.zeros(tr.shape, N.float_)
            t    = tr[:, 0:1]
            r    = tr[:, 1:2]
            x    = xy[:, 0:1]
            y    = xy[:, 1:2]
            x[:] = r * N.sin(t)
            y[:] = r * N.cos(t)
            return xy

        transform_non_affine = transform

        def inverted(self):
            return NorthPolarAxes.InvertedNorthPolarTransform()

    class InvertedNorthPolarTransform(MPLproj.PolarAxes.InvertedPolarTransform):
        def transform(self, xy):
            x = xy[:, 0:1]
            y = xy[:, 1:]
            r = N.sqrt(x*x + y*y)
            theta = N.arctan2(y, x)
            return N.concatenate((theta, r), 1)

        def inverted(self):
            return NorthPolarAxes.NorthPolarTransform()

    def _set_lim_and_transforms(self):
        MPLproj.PolarAxes._set_lim_and_transforms(self)
        self.transProjection = self.NorthPolarTransform()
        self.transData = (
            self.transScale +
            self.transProjection +
            (self.transProjectionAffine + self.transAxes))
        self._xaxis_transform = (
            self.transProjection +
            self.PolarAffine(MPLtrans.IdentityTransform(),
                             MPLtrans.Bbox.unit()) + self.transAxes)
        self._xaxis_text1_transform = (
            self._theta_label1_position +
            self._xaxis_transform)
        self._yaxis_transform = (
            MPLtrans.Affine2D().scale(N.pi * 2.0, 1.0) +
            self.transData)
        self._yaxis_text1_transform = (
            self._r_label1_position +
            MPLtrans.Affine2D().scale(1.0 / 360.0, 1.0) +
            self._yaxis_transform)

if __name__ == "__main__":
  register_projection(NorthPolarAxes)

  angle = N.arange(0, 360, 10, dtype=float) * N.pi / 180.0
  arbitrary_data = (N.abs(N.sin(angle)) + 0.1 *
    (N.random.random_sample(size=angle.shape) - 0.5))

  P.clf()
  P.subplot(1, 1, 1, projection='northpolar')
  P.plot(angle, arbitrary_data)
  P.show()
 
