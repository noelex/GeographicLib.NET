using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace GeographicLib
{
    static class Extensions
    {
        public static bool LessThan<T>(this T t, T other) where T : IComparable<T> => t.CompareTo(other) < 0;

        public static bool LessThanOrEqualTo<T>(this T t, T other) where T : IComparable<T> => t.CompareTo(other) <= 0;

        public static bool EqualTo<T>(this T t, T other) where T : IComparable<T> => t.CompareTo(other) == 0;

        public static bool GreaterThan<T>(this T t, T other) where T : IComparable<T> => t.CompareTo(other) > 0;

        public static bool GreaterThanOrEqualTo<T>(this T t, T other) where T : IComparable<T> => t.CompareTo(other) >= 0;

        public static void Swap<T>(this IList<T> list, int indexA, int indexB)
        {
            var t = list[indexA];
            list[indexA] = list[indexB];
            list[indexB] = t;
        }

        public static void NthElement<T>(this IList<T> array, int startIndex, int nthSmallest, int endIndex, IComparer<T> comparer)
        {
            int from = startIndex;
            int to = endIndex;

            // if from == to we reached the kth element
            while (from < to)
            {
                int r = from, w = to;
                T mid = array[(r + w) / 2];

                // stop if the reader and writer meets
                while (r < w)
                {
                    if (comparer.Compare(array[r], mid) > -1)
                    { // put the large values at the end
                        T tmp = array[w];
                        array[w] = array[r];
                        array[r] = tmp;
                        w--;
                    }
                    else
                    { // the value is smaller than the pivot, skip
                        r++;
                    }
                }

                // if we stepped up (r++) we need to step one down
                if (comparer.Compare(array[r], mid) > 0)
                {
                    r--;
                }

                // the r pointer is on the end of the first k elements
                if (nthSmallest <= r)
                {
                    to = r;
                }
                else
                {
                    from = r + 1;
                }
            }

            return;
        }
    }

    internal class NearestNeighbor<TPoint>
    {
        // For tracking changes to the I/O format.
        private const int version = 1;

        // This is what we get "free"; but if sizeof(dist_t) = 1 (unlikely), allow
        // 4 slots (and this accommodates the default value bucket = 4).
        private const int maxbucket = (2 + ((4 * sizeof(double)) / sizeof(int) >= 2 ?
            (4 * sizeof(double)) / sizeof(int) : 2));

        private static readonly ItemComparer Comparer = new ItemComparer();

        private int _numpoints, _bucket, _cost;
        private List<Node> _tree;

        // Counters to track stastistics on the cost of searches
        private double _mc, _sc;
        private int _c1, _k, _cmin, _cmax;

        public NearestNeighbor()
        {
            (_numpoints, _bucket, _cost) = (0, 0, 0);
        }

        public NearestNeighbor(IList<TPoint> pts, Func<TPoint, TPoint, double> distfunct, int bucket = 4)
        {
            Initialize(pts, distfunct, bucket);
        }

        public int NumPoints => _numpoints;

        public void Initialize(IList<TPoint> pts, Func<TPoint, TPoint, double> dist,
                    int bucket = 4)
        {
            if (!(0 <= bucket && bucket <= maxbucket))
                throw new GeographicException("bucket must lie in [0, 2 + 4*sizeof(dist_t)/sizeof(int)]");

            // the pair contains distance+id
            var ids = Enumerable.Range(0, pts.Count())
                .Select(k => new item(0,k)).ToList();

            int cost = 0;
            _tree = new List<Node>();
            Init(pts, dist, bucket, _tree, ids, ref cost,
                 0, ids.Count, ids.Count / 2);

            _numpoints = pts.Count;
            _bucket = bucket;
            _mc = _sc = 0;
            _cost = cost; _c1 = _k = _cmax = 0;
            _cmin = int.MaxValue;
        }

        public double Search(IList<TPoint> pts, Func<TPoint, TPoint, double> dist, TPoint query,
            List<int> ind, int k = 1, double maxdist = double.MaxValue, double mindist = -1, bool exhaustive = true, double tol = 0)
        {
            if (_numpoints != pts.Count)
                throw new GeographicException("pts array has wrong size");

            double d;
            var results = new PriorityQueue<item>(Comparer, false);
            if (_numpoints > 0 && k > 0 && maxdist > mindist)
            {
                // distance to the kth closest point so far
                var tau = maxdist;
                // first is negative of how far query is outside boundary of node
                // +1 if on boundary or inside
                // second is node index
                var todo = new PriorityQueue<item>(Comparer, false);
                todo.Enqueue(new item(1, _tree.Count - 1));
                int c = 0;

                while (todo.Any())
                {
                    int n = todo.First().second;
                    d = -todo.First().first;
                    todo.Dequeue();
                    var tau1 = tau - tol;
                    // compare tau and d again since tau may have become smaller.
                    if (!(n >= 0 && tau1 >= d)) continue;
                    var current = _tree[n];
                    var dst = 0.0;   // to suppress warning about uninitialized variable
                    bool exitflag = false, leaf = current.Index < 0;
                    for (int i = 0; i < (leaf ? _bucket : 1); ++i)
                    {
                        int index = leaf ? current.Leaves.Span[i] : current.Index;
                        if (index < 0) break;
                        dst = dist(pts[index], query);
                        ++c;

                        if (dst > mindist && dst <= tau)
                        {
                            if (results.Count == k) results.Dequeue();
                            results.Enqueue(new item(dst, index));
                            if (results.Count == k)
                            {
                                if (exhaustive)
                                    tau = results.First().first;
                                else
                                {
                                    exitflag = true;
                                    break;
                                }
                                if (tau <= tol)
                                {
                                    exitflag = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (exitflag) break;

                    if (current.Index < 0) continue;
                    tau1 = tau - tol;
                    if (current.Data.child0 >= 0 && dst + current.Data.upper0 >= mindist)
                    {
                        if (dst < current.Data.lower0)
                        {
                            d = current.Data.lower0 - dst;
                            if (tau1 >= d)
                                todo.Enqueue(new item(-d, current.Data.child0));
                        }
                        else if (dst > current.Data.upper0)
                        {
                            d = dst - current.Data.upper0;
                            if (tau1 >= d)
                                todo.Enqueue(new item(-d, current.Data.child0));
                        }
                        else
                            todo.Enqueue(new item(1, current.Data.child0));
                    }

                    if (current.Data.child1 >= 0 && dst + current.Data.upper1 >= mindist)
                    {
                        if (dst < current.Data.lower1)
                        {
                            d = current.Data.lower1 - dst;
                            if (tau1 >= d)
                                todo.Enqueue(new item(-d, current.Data.child1));
                        }
                        else if (dst > current.Data.upper1)
                        {
                            d = dst - current.Data.upper1;
                            if (tau1 >= d)
                                todo.Enqueue(new item(-d, current.Data.child1));
                        }
                        else
                            todo.Enqueue(new item(1, current.Data.child1));
                    }
                }
                ++_k;
                _c1 += c;
                double omc = _mc;
                _mc += (c - omc) / _k;
                _sc += (c - omc) * (c - _mc);
                if (c > _cmax) _cmax = c;
                if (c < _cmin) _cmin = c;
            }

            d = -1.0;

            for (var i = 0; i < results.Count; i++)
            {
                ind.Add(results.First().second);
                if (i == 0) d = results.First().first;
                results.Dequeue();
            }

            return d;
        }

        public SearchStatistics Statistics =>
            new SearchStatistics
            { 
                SetupCost = _cost,
                NumSearches = _k,
                SearchCost = _c1,
                MinCost = _cmin,
                MaxCost = _cmax,
                Mean = _mc,
                Sd = Math.Sqrt(_sc / (_k - 1))
            };

        public void ResetStatistics()
        {
            _mc = _sc = 0;
            _c1 = _k = _cmax = 0;
            _cmin = int.MaxValue;
        }

        private int Init(IList<TPoint> pts, Func<TPoint, TPoint, double> dist, int bucket,
                IList<Node> tree, List<item> ids, ref int cost, int l, int u, int vp)
        {
            if (u == l)
                return -1;
            var node = new Node();

            if (u - l > (bucket == 0 ? 1 : bucket))
            {

                // choose a vantage point and move it to the start
                int i = vp;
                ids.Swap(l, i);

                int m = (u + l + 1) / 2;

                for (int k = l + 1; k < u; ++k)
                {
                    ids[k].first = dist(pts[ids[l].second], pts[ids[k].second]);
                    ++cost;
                }

                ids.NthElement(l + 1, m, u, Comparer);

                node.Index = ids[l].second;

                item t;
                int ti;

                if (m > l + 1)
                {        // node.child[0] is possibly empty
                    t = ids.Skip(l + 1).Take(m - (l + 1)).OrderBy(x => x.first).First();
                    node.Data.lower0 = t.first;
                    t = ids.Skip(l + 1).Take(m - (l + 1)).OrderByDescending(x => x.first).First();
                    ti = ids.IndexOf(t);
                    node.Data.upper0 = t.first;
                    // Use point with max distance as vantage point; this point act as a
                    // "corner" point and leads to a good partition.
                    node.Data.child0 = Init(pts, dist, bucket, tree, ids, ref cost,
                                              l + 1, m, ti);
                }

                t = ids.Skip(m).Take(u - m).OrderByDescending(x => x.first).First();
                ti = ids.IndexOf(t);
                node.Data.lower1 = ids[m].first;
                node.Data.upper1 = t.first;
                // Use point with max distance as vantage point here too
                node.Data.child1 = Init(pts, dist, bucket, tree, ids, ref cost,
                                          m, u, ti);
            }
            else
            {
                if (bucket == 0)
                    node.Index = ids[l].second;
                else
                {
                    node.Index = -1;
                    // Sort the bucket entries so that the tree is independent of the
                    // implementation of nth_element.
                    ids.Sort(l, u - l, Comparer);

                    for (int i = l; i < u; ++i)
                        node.Leaves.Span[i - l] = ids[i].second;
                    for (int i = u - l; i < bucket; ++i)
                        node.Leaves.Span[i] = -1;
                    for (int i = bucket; i < maxbucket; ++i)
                        node.Leaves.Span[i] = 0;
                }
            }

            tree.Add(node);
            return tree.Count - 1;
        }

        private class Node
        {
            private Bounds _data;

            public Node()
            {
                _data = MemoryMarshal.Cast<int, Bounds>(Leaves.Span)[0];
                _data.child0 = _data.child1 = -1;
            }

            public Memory<int> Leaves { get; } = new int[maxbucket];

            public ref Bounds Data => ref _data;

            public int Index { get; set; } = -1;

            public void Check(int numpoints, int treesize, int bucket)
            {
                if (!(-1 <= Index && Index < numpoints))
                    throw new GeographicException("Bad index");
                if (Index >= 0)
                {
                    if (!(-1 <= _data.child0 && _data.child0 < treesize &&
                           -1 <= _data.child1 && _data.child1 < treesize))
                        throw new GeographicException("Bad child pointers");
                    if (!(0 <= _data.lower0 && _data.lower0 <= _data.upper0 &&
                           _data.upper0 <= _data.lower1 &&
                           _data.lower1 <= _data.upper1))
                        throw new GeographicException("Bad bounds");
                }
                else
                {
                    // Must be at least one valid leaf followed by a sequence end markers
                    // (-1).
                    bool start = true;
                    for (int l = 0; l < bucket; ++l)
                    {
                        if (!((start ?
                                ((l == 0 ? 0 : -1) <= Leaves.Span[l] && Leaves.Span[l] < numpoints) :
                                Leaves.Span[l] == -1)))
                            throw new GeographicException("Bad leaf data");
                        start = Leaves.Span[l] >= 0;
                    }
                    for (int l = bucket; l < maxbucket; ++l)
                    {
                        if (Leaves.Span[l] != 0)
                            throw new GeographicException("Bad leaf data");
                    }
                }
            }

            public struct Bounds
            {
                public double lower0, lower1, upper0, upper1; // bounds on inner/outer distances
                public int child0, child1;
            };
        }

        private class item
        {
            public double first;
            public int second;

            public item(double f, int s)
            {
                first = f;
                second = s;
            }
        }

        private class ItemComparer : IComparer<item>
        {
            public int Compare(item x, item y)
            {
                return x.first.CompareTo(y.first);
            }
        }

        public class SearchStatistics
        {
            public int SetupCost { get; internal set; }

            public int NumSearches { get; internal set; }

            public int SearchCost { get; internal set; }

            public int MinCost { get; internal set; }

            public int MaxCost { get; internal set; }

            public double Mean { get; internal set; }

            public double Sd { get; internal set; }
        }
    }
}
