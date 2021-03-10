using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    public class PolygonArea<T> where T : IGeodesicLike
    {
        private readonly T _earth;

        private double _area0;
        private bool _polyline;
        private uint _mask, _num;
        private int _crossings;

    }
}
