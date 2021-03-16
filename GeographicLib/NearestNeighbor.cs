using System;
using System.Collections.Generic;
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

        public static void Swap<T>(this List<T> list, int indexA,int indexB)
        {
            var t = list[indexA];
            list[indexA] = list[indexB];
            list[indexB] = t;
        }

        public static void NthElement<T>(this List<T> array, int startIndex, int nthSmallest, int endIndex, IComparer<T> comparer)
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

    public class NearestNeighbor<TDistance, TPoint>
        where TDistance : struct, IComparable<TDistance>
    {
        // For tracking changes to the I/O format.
        private const int version = 1;

        private static readonly TDistance Zero = default(TDistance);

        private static readonly ItemComparer Comparer = new ItemComparer();

        // This is what we get "free"; but if sizeof(dist_t) = 1 (unlikely), allow
        // 4 slots (and this accommodates the default value bucket = 4).
        private static readonly int maxbucket = (2 + ((4 * Marshal.SizeOf<TDistance>()) / sizeof(int) >= 2 ?
            (4 * Marshal.SizeOf<TDistance>()) / sizeof(int) : 2));

        private readonly int _numpoints, _bucket, _cost;
        private readonly List<Node> _tree = new List<Node>();

        // Counters to track stastistics on the cost of searches
        private double _mc, _sc;
        private int _c1, _k, _cmin, _cmax;

        private int Init(List<TPoint> pts, Func<TPoint,TPoint,TDistance> dist, int bucket,
            List<Node> tree, List<item> ids, ref int cost, int l, int u, int vp)
        {
            if (u == l)
                return -1;
            var node = new Node();

            if (u - l > (bucket == 0 ? 1 : bucket))
            {

                // choose a vantage point and move it to the start
                int i = vp;
                ids.Swap(l,i);

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
                node.Data.child1 = Init(pts, dist, bucket, tree, ids,ref cost,
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
                    if (!(Zero.LessThanOrEqualTo(_data.lower0) && _data.lower0.LessThanOrEqualTo(_data.upper0) &&
                           _data.upper0.LessThanOrEqualTo(_data.lower1) &&
                           _data.lower1.LessThanOrEqualTo(_data.upper1)))
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
                public TDistance lower0, lower1, upper0, upper1; // bounds on inner/outer distances
                public int child0, child1;
            };
        }

        private class item
        {
            public TDistance first;
            public int second;
        }

        private class ItemComparer : IComparer<item>
        {
            public int Compare(item x, item y)
            {
                return x.first.CompareTo(y.first);
            }
        }
    }
}
