pub mod counting;
pub mod proto;

use crate::bloomfilter::proto::Proto;
use bitvector::*;
use std::rc::Rc;
use std::vec::Vec;

type BloomFilterType = Box<dyn BloomFilter>;

/// The interface between the internal and external representations of the contents of a Bloom filter.

/// The traits that all BloomFilters must share.
///
/// # Note
/// All Bloom filters know their shape.  However, default function implementations ignore shape
/// differences before executing comparisons and other multi-filter operations.  This is because
/// in most cases all filters will have the same shape.  In cases where it is important to check
/// it is up to the developer to check the `shape.equivalent_to()` method before making the calls.
///
pub trait BloomFilter {
    /// Merges two Bloom filters to create a new one.  In the default implementation, if the
    /// merge will probably produce a sparse filter the sparse
    /// implementation is used.  Otherwise the Simple implementation is used.
    /// If the bloom filter arguments are not of Simple or Sparse they will be treated as though
    /// they were.  (e.g. no special handling for other types)
    ///
    /// Other implementations may exist in other BloomFilter implementations.  Those implementations
    /// are expected to produce new instances of their type.
    ///
    /// # Arguments
    /// * `other` - The other bloom filter to merge
    fn merge(&self, other: &BloomFilterType) -> Result<BloomFilterType, &str> {
        if (self.hamming_value() + other.hamming_value()) < self.shape().number_of_buckets() {
            // create a sparse one
            print!("Creating sparse BloomFilter\n");
            let mut x = Vec::with_capacity(self.hamming_value() + other.hamming_value());
            x.extend(self.indices().iter());
            x.extend(other.indices().iter());
            x.sort_unstable();
            x.dedup();
            Ok(Box::new(Sparse {
                shape: self.shape().clone(),
                buffer: Rc::new(x),
            }))
        } else {
            // create a simple one
            print!("Creating simple BloomFilter\n");
            Ok(Box::new(Simple {
                shape: self.shape().clone(),
                buffer: Rc::new(self.vector().union(&other.vector())),
            }))
        }
    }

    /// Updates this filter by merging the values from the other filter.
    /// This is an atomic function in that it either works or it fails, but if it fails
    /// the internal structure is as it was before the method was called.
    ///
    /// # Arguments
    /// * `other` - The other bloom filter to merge into this one.
    fn merge_inplace<'a>(&mut self, other: &BloomFilterType) -> Result<(), &'a str>;

    /// Merges a proto Bloom filter into a bloom filter.  In the default implementation, if the
    /// merge will probably produce a sparse filter the sparse
    /// implementation is used.  Otherwise the Simple implementation is used.
    ///
    /// Other implementations may exist in other BloomFilter implementations.  Those implementations
    /// are expected to produce new instances of their type.
    ///
    /// # Arguments
    /// * `proto` - The proto bloom filter to merge into this one.
    fn merge_proto(&self, proto: &dyn Proto) -> Result<BloomFilterType, &str> {
        let other = BloomFilterFactory::materialize(self.shape(), proto);
        return self.merge(&other);
    }

    /// Merges a proto Bloom filter into this bloom filter.
    /// This is an atomic function in that it either works or it fails, but if it fails
    /// the internal structure is as it was before the method was called.
    ///
    /// # Arguments
    /// * `proto` - The proto bloom filter to merge into this one.
    fn merge_proto_inplace<'a>(&mut self, proto: &dyn Proto) -> Result<(), &'a str> {
        let other = BloomFilterFactory::materialize(self.shape(), proto);
        return self.merge_inplace(&other);
    }

    /// Returns true if it is more efficient to access the bits by index rather than by BitVector.
    /// In the general the case this is true if hamming_value < shape.m.  For other implementations
    /// it may be easier to produce a representation as a list of index values (sparse) or a
    /// BitVector (not sparse).
    fn is_sparse(&self) -> bool {
        self.shape().is_sparse(self.hamming_value())
    }

    /// Returns a the BitVector that represents this bloom filter.
    fn vector(&self) -> Rc<BitVector>;

    /// Returns a Vector of the bits that are turned on.
    /// these are required to be in sorted into ascending order and duplicates are not allowed.
    fn indices(&self) -> Rc<Vec<usize>>;

    /// Returns the shape of the filter
    fn shape(&self) -> &Shape;

    /// Gets the hamming value (the number of bits turned on).
    fn hamming_value(&self) -> usize;

    /// Determines if this filter contains the proto Bloom filter.
    ///
    /// # Arguments
    /// * `proto` - The proto bloom filter to check.
    fn contains_proto(&self, proto: &dyn Proto) -> bool {
        let other = BloomFilterFactory::materialize(self.shape(), proto);
        return self.contains(&other);
    }

    /// Determines if this filter contains the other filter.
    ///
    /// # Arguments
    /// * `other` - The other bloom filter to check.
    fn contains(&self, other: &BloomFilterType) -> bool {
        // this is counter intuitive.  We skip all the matches and find the first
        // non matching BitBucket.  If after this there next() return None then the
        // the filters match.
        return if self.is_sparse() {
            print!("self is sparse, ");
            if other.is_sparse() {
                print!(" other is sparse\n");
                if other.indices().len() > self.indices().len() {
                    false
                } else {
                    other
                        .indices()
                        .iter()
                        .all(|s| self.indices().binary_search(&s).is_ok())
                }
            } else {
                print!("other is not sparse\n");
                let x: &BitVector = &other.vector();
                x.into_iter()
                    .all(|s| self.indices().binary_search(&s).is_ok())
            }
        } else {
            print!("self is not sparse, ");
            if other.is_sparse() {
                print!("other is sparse\n");
                other.indices().iter().all(|s| self.vector().contains(*s))
            } else {
                print!("other is not sparse\n");
                let v1: &BitVector = &self.vector();
                let v2: &BitVector = &other.vector();
                (v1 & v2).len() == v2.len()
            }
        };
    }

    //
    //  SIZE METHODS
    //

    /// Estimates the number of items in this filter.
    fn estimate_n(&self) -> f64 {
        self.shape().estimate_n(self.hamming_value())
    }

    /// Estimates the number of items in the union of the two filters.
    ///
    /// # Arguments
    /// * `other` - The other bloom filter in the union.
    fn estimate_union(&self, other: &BloomFilterType) -> f64 {
        self.shape()
            .estimate_n(self.vector().union(&other.vector()).len())
    }

    /// Estimates the number of items in the intersection of the two filters.
    ///
    /// # Arguments
    /// * `other` - The other bloom filter in the intersection.
    fn estimate_intersection(&self, other: &BloomFilterType) -> f64 {
        self.estimate_n() + other.estimate_n() - self.estimate_union(other)
    }
}

/// The shape of the bloom filter.
#[derive(Copy, Clone, Debug)]
pub struct Shape {
    /// number of functions
    k: usize,
    /// number of bits in the filter
    m: usize,
}

impl Shape {
    ///  Gets the number of u64 sized buckets necessary for to represent the bits.
    pub fn number_of_buckets(&self) -> usize {
        let mut len = self.m / 64;
        if self.m % 64 > 0 {
            len += 1;
        }
        len
    }

    /// Determines if a bloom filter with hamming bits enabled would be considered "sparse" or not.
    /// A sparse filter is one in which there are likely to be empty buckets.
    ///
    /// # Arguments
    /// * `hamming` - The hamming value to check.
    pub fn is_sparse(&self, hamming: usize) -> bool {
        hamming < self.number_of_buckets()
    }

    /// Checks if another shape is equivalent to this one.  Shapes are equivalent if the `m` and `k`
    /// values are the same.
    ///
    /// # Arguments
    /// * `other` - The other Shape to check against.
    pub fn equivalent_to(&self, other: &Shape) -> bool {
        self.m != other.m || self.k != other.k
    }

    /// Calculates the probability of false positives if `count` objects are added to a filter
    /// using this shape.
    ///
    /// # Arguments
    /// * `count` - The number of items in the Bloom filter.
    pub fn false_positives(&self, count: usize) -> f64 {
        let k = self.k as f64;
        let n = count as f64;
        let m = self.m as f64;
        (1.0 - ((-k * n) / m).exp()).powf(k)
    }

    /// Calculates the estimated number of items in a filter
    /// of this Shape.
    ///
    /// # Arguments
    /// * `count` - - he number of bits turned on in the filter.
    ///
    pub fn estimate_n(&self, count: usize) -> f64 {
        let c = count as f64;
        let m = self.m as f64;
        let k = self.k as f64;
        let result = -(m / k) * (1.0 - (c / m)).ln();
        result
    }
}

/// A bloom filter that stores the bits in a BitVector
#[derive(Debug)]
pub struct Simple {
    shape: Shape,
    buffer: Rc<BitVector>,
}

impl Simple {
    /// Create an empty Simple instance
    ///
    /// # Arguments
    /// * `shape` - The shape of the Bloom filter.
    pub fn empty_instance(shape: &Shape) -> BloomFilterType {
        Box::new(Simple {
            shape: shape.clone(),
            buffer: Rc::new(BitVector::new(shape.number_of_buckets())),
        })
    }

    /// Create a Simple instance from a proto type.
    ///
    /// # Arguments
    /// * `shape` - The shape of the Bloom filter.
    /// * `proto` - The prototype for the filter.
    pub fn instance(shape: &Shape, proto: &dyn Proto) -> BloomFilterType {
        let mut v = BitVector::new(shape.number_of_buckets());
        proto.bits(&shape).into_iter().for_each(|s| {
            v.insert(s);
        });
        let simple = Simple {
            shape: shape.clone(),
            buffer: Rc::new(v),
        };
        return Box::new(simple);
    }
}

impl BloomFilter for Simple {
    fn merge_inplace<'a>(&mut self, other: &BloomFilterType) -> Result<(), &'a str> {
        let mut b = BitVector::new(self.shape.number_of_buckets());
        b.insert_all(&self.buffer);
        if other.is_sparse() {
            other.indices().iter().for_each(|s| {
                b.insert(*s);
            });
        } else {
            b.insert_all(&other.vector());
        }
        self.buffer = Rc::new(b);
        Ok(())
    }

    /// Always returns `false`
    fn is_sparse(&self) -> bool {
        false
    }

    /// return the filter as a BitVector.
    fn vector(&self) -> Rc<BitVector> {
        self.buffer.clone()
    }

    /// return a list of the bits that are turned on.
    fn indices(&self) -> Rc<Vec<usize>> {
        let x: &BitVector = &self.buffer;
        Rc::new(x.into_iter().collect())
    }

    fn shape(&self) -> &Shape {
        return &self.shape;
    }

    fn hamming_value(&self) -> usize {
        return self.buffer.len();
    }
}

/// A bloom filter that stores the index of the enabled bits.  Useful for Bloomfilters of shapes
/// with large `m` and small `k` values.  Particularly where `k` < m/64.
///
/// # See
/// * Shape.is_sparse()
#[derive(Debug)]
pub struct Sparse {
    shape: Shape,
    buffer: Rc<Vec<usize>>,
}

impl Sparse {
    /// Create an empty Sparse instance.
    ///
    /// # Arguments
    /// * `shape` - The shape of the Bloom filter.
    pub fn empty_instance(shape: &Shape) -> BloomFilterType {
        Box::new(Sparse {
            shape: shape.clone(),
            buffer: Rc::new(Vec::new()),
        })
    }

    /// Create a Sparse instance from a proto type.
    ///
    /// # Arguments
    /// * `shape` - The shape of the Bloom filter.
    /// * `proto` - The prototype for the filter.
    pub fn instance(shape: &Shape, proto: &dyn Proto) -> BloomFilterType {
        let mut v: Vec<usize> = proto.bits(&shape).into_iter().collect();
        v.sort_unstable();
        v.dedup();
        let sparse = Sparse {
            shape: shape.clone(),
            buffer: Rc::new(v),
        };
        return Box::new(sparse);
    }
}

impl BloomFilter for Sparse {
    fn merge_inplace<'a>(&mut self, other: &BloomFilterType) -> Result<(), &'a str> {
        let mut v: Vec<usize> = self.buffer.to_vec();
        if other.is_sparse() {
            other.indices().iter().for_each(|s| v.push(*s));
        } else {
            let x: &BitVector = &other.vector();
            x.into_iter().for_each(|s| v.push(s));
        }
        v.sort_unstable();
        v.dedup();
        self.buffer = Rc::new(v);
        Ok(())
    }

    /// Always returns `true`
    fn is_sparse(&self) -> bool {
        true
    }

    /// return the filter as a BitVector.
    fn vector(&self) -> Rc<BitVector> {
        let mut v = BitVector::new(self.shape.number_of_buckets());
        self.buffer.iter().for_each(|s| {
            v.insert(*s);
        });
        Rc::new(v)
    }

    /// return a list of the bits that are turned on.
    fn indices(&self) -> Rc<Vec<usize>> {
        self.buffer.clone()
    }

    fn shape(&self) -> &Shape {
        return &self.shape;
    }

    fn hamming_value(&self) -> usize {
        return self.buffer.len();
    }
}

#[cfg(test)]
mod tests {

    use crate::bloomfilter::proto::*;
    use crate::bloomfilter::*;

    #[test]
    fn shape_false_positives() {
        let shape = Shape { m: 134_191, k: 23 };
        assert!(shape.false_positives(4000) - (1.0 / 9_994_297.0) < 0.0000001);
    }

    #[test]
    fn shape_number_of_buckets() {
        let shape = Shape { m: 60, k: 2 };
        assert_eq!(shape.number_of_buckets(), 1);
        let shape = Shape { m: 120, k: 2 };
        assert_eq!(shape.number_of_buckets(), 2);
    }

    #[test]
    fn empty_filter() {
        let shape = Shape { m: 60, k: 2 };
        let bloomfilter = Simple::empty_instance(&shape);
        assert_eq!(*bloomfilter.indices(), []);
        let v = bloomfilter.vector();
        assert_eq!(v.len(), 0);
        assert_eq!(bloomfilter.hamming_value(), 0);
        assert!(bloomfilter.estimate_n() < 0.05);
        // filter always contains itself
        assert!(bloomfilter.contains(&bloomfilter));
    }

    #[test]
    fn filter_build_correct() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let bloomfilter = Simple::instance(&shape, &proto);
        assert_eq!(*bloomfilter.indices(), [0, 1]);
        let v = bloomfilter.vector();
        assert_eq!(v.len(), 2);
        assert!(v.contains(0));
        assert!(v.contains(1));
        assert_eq!(bloomfilter.hamming_value(), 2);
        assert!(bloomfilter.estimate_n() - 1.0 < 0.05);
        // filter always contains itself
        assert!(bloomfilter.contains(&bloomfilter));
        let empty_filter: BloomFilterType = Box::new(Simple {
            shape: shape.clone(),
            buffer: Rc::new(BitVector::new(shape.number_of_buckets())),
        });
        // a filter always contains the empty filter
        assert!(bloomfilter.contains(&empty_filter));
        // an empty filter never contains a populated filter
        print!("{:#?}", empty_filter.contains(&bloomfilter));
        assert!(!empty_filter.contains(&bloomfilter));
        // an empty filter always contains itself
        assert!(empty_filter.contains(&empty_filter));
    }

    #[test]
    fn contains_test_simple_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let proto2 = SimpleProto::from_u32(5);
        let bloomfilter = Simple::instance(&shape, &proto);

        let mut bloomfilter2 = Simple::instance(&shape, &proto);
        assert!(bloomfilter2
            .merge_inplace(&Simple::instance(&shape, &proto2))
            .is_ok());

        let simple_empty_filter: BloomFilterType = Box::new(Simple {
            shape: shape.clone(),
            buffer: Rc::new(BitVector::new(shape.number_of_buckets())),
        });

        let sparse_empty_filter: BloomFilterType = Box::new(Sparse {
            shape: shape.clone(),
            buffer: Rc::new(vec![]),
        });

        // filter always contains itself
        assert!(bloomfilter.contains(&bloomfilter));
        // a filter always contains the empty filter
        assert!(bloomfilter.contains(&simple_empty_filter));
        assert!(bloomfilter.contains(&sparse_empty_filter));
        // an empty filter never contains a populated filter
        assert!(!simple_empty_filter.contains(&bloomfilter));
        assert!(!sparse_empty_filter.contains(&bloomfilter));
        // an empty filter always contains itself
        assert!(simple_empty_filter.contains(&simple_empty_filter));
        assert!(simple_empty_filter.contains(&sparse_empty_filter));
        assert!(sparse_empty_filter.contains(&simple_empty_filter));
        assert!(sparse_empty_filter.contains(&sparse_empty_filter));

        // bloom filter2 contains bloom filter 1
        assert!(bloomfilter2.contains(&bloomfilter));
        // bloom filter1 does not contain bloom filter 2
        assert!(!bloomfilter.contains(&bloomfilter2));
    }

    #[test]
    fn contains_test_simple_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let proto2 = SimpleProto::from_u32(5);
        let bloomfilter = Simple::instance(&shape, &proto);

        let mut bloomfilter2 = Sparse::instance(&shape, &proto);
        assert!(bloomfilter2
            .merge_inplace(&Simple::instance(&shape, &proto2))
            .is_ok());

        let simple_empty_filter: BloomFilterType = Box::new(Simple {
            shape: shape.clone(),
            buffer: Rc::new(BitVector::new(shape.number_of_buckets())),
        });

        let sparse_empty_filter: BloomFilterType = Box::new(Sparse {
            shape: shape.clone(),
            buffer: Rc::new(vec![]),
        });

        // filter always contains itself
        assert!(bloomfilter.contains(&bloomfilter));
        // a filter always contains the empty filter
        assert!(bloomfilter.contains(&simple_empty_filter));
        assert!(bloomfilter.contains(&sparse_empty_filter));
        // an empty filter never contains a populated filter
        assert!(!simple_empty_filter.contains(&bloomfilter));
        assert!(!sparse_empty_filter.contains(&bloomfilter));
        // an empty filter always contains itself
        assert!(simple_empty_filter.contains(&simple_empty_filter));
        assert!(simple_empty_filter.contains(&sparse_empty_filter));
        assert!(sparse_empty_filter.contains(&simple_empty_filter));
        assert!(sparse_empty_filter.contains(&sparse_empty_filter));

        // bloom filter2 contains bloom filter 1
        assert!(bloomfilter2.contains(&bloomfilter));
        // bloom filter1 does not contain bloom filter 2
        assert!(!bloomfilter.contains(&bloomfilter2));
    }

    #[test]
    fn contains_test_sparse_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let proto2 = SimpleProto::from_u32(5);
        let bloomfilter = Sparse::instance(&shape, &proto);

        let mut bloomfilter2 = Simple::instance(&shape, &proto);
        assert!(bloomfilter2
            .merge_inplace(&Simple::instance(&shape, &proto2))
            .is_ok());

        let simple_empty_filter: BloomFilterType = Box::new(Simple {
            shape: shape.clone(),
            buffer: Rc::new(BitVector::new(shape.number_of_buckets())),
        });

        let sparse_empty_filter: BloomFilterType = Box::new(Sparse {
            shape: shape.clone(),
            buffer: Rc::new(vec![]),
        });

        // filter always contains itself
        assert!(bloomfilter.contains(&bloomfilter));
        // a filter always contains the empty filter
        assert!(bloomfilter.contains(&simple_empty_filter));
        assert!(bloomfilter.contains(&sparse_empty_filter));
        // an empty filter never contains a populated filter
        assert!(!simple_empty_filter.contains(&bloomfilter));
        assert!(!sparse_empty_filter.contains(&bloomfilter));
        // an empty filter always contains itself
        assert!(simple_empty_filter.contains(&simple_empty_filter));
        assert!(simple_empty_filter.contains(&sparse_empty_filter));
        assert!(sparse_empty_filter.contains(&simple_empty_filter));
        assert!(sparse_empty_filter.contains(&sparse_empty_filter));

        // bloom filter2 contains bloom filter 1
        assert!(bloomfilter2.contains(&bloomfilter));
        // bloom filter1 does not contain bloom filter 2
        assert!(!bloomfilter.contains(&bloomfilter2));
    }

    #[test]
    fn contains_test_sparse_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let proto2 = SimpleProto::from_u32(5);
        let bloomfilter = Sparse::instance(&shape, &proto);

        let mut bloomfilter2 = Sparse::instance(&shape, &proto);
        assert!(bloomfilter2
            .merge_inplace(&Simple::instance(&shape, &proto2))
            .is_ok());

        let simple_empty_filter: BloomFilterType = Box::new(Simple {
            shape: shape.clone(),
            buffer: Rc::new(BitVector::new(shape.number_of_buckets())),
        });

        let sparse_empty_filter: BloomFilterType = Box::new(Sparse {
            shape: shape.clone(),
            buffer: Rc::new(vec![]),
        });

        // filter always contains itself
        assert!(bloomfilter.contains(&bloomfilter));
        // a filter always contains the empty filter
        assert!(bloomfilter.contains(&simple_empty_filter));
        assert!(bloomfilter.contains(&sparse_empty_filter));
        // an empty filter never contains a populated filter
        assert!(!simple_empty_filter.contains(&bloomfilter));
        assert!(!sparse_empty_filter.contains(&bloomfilter));
        // an empty filter always contains itself
        assert!(simple_empty_filter.contains(&simple_empty_filter));
        assert!(simple_empty_filter.contains(&sparse_empty_filter));
        assert!(sparse_empty_filter.contains(&simple_empty_filter));
        assert!(sparse_empty_filter.contains(&sparse_empty_filter));

        // bloom filter2 contains bloom filter 1
        assert!(bloomfilter2.contains(&bloomfilter));
        // bloom filter1 does not contain bloom filter 2
        assert!(!bloomfilter.contains(&bloomfilter2));
    }

    #[test]
    fn shape_used_multiple_times() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let bloomfilter = Simple::instance(&shape, &proto);
        let bloomfilter2 = Simple::instance(&shape, &proto);
        assert_eq!(*bloomfilter.indices(), [0, 1]);
        assert_eq!(*bloomfilter2.indices(), [0, 1]);
        assert!(bloomfilter.contains(&bloomfilter2));
    }

    #[test]
    fn filter_merge_inplace_test_simple_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let mut bloomfilter = Simple::instance(&shape, &proto);
        let proto2 = SimpleProto::from_u32(0x100);
        let bloomfilter2 = Simple::instance(&shape, &proto2);

        assert!(bloomfilter.merge_inplace(&bloomfilter2).is_ok());
        assert_eq!(*bloomfilter.indices(), [0, 1, 16]);

        assert!(bloomfilter.contains(&bloomfilter2));
    }

    #[test]
    fn filter_merge_inplace_test_simple_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let mut bloomfilter = Simple::instance(&shape, &proto);
        let proto2 = SimpleProto::from_u32(0x100);
        let bloomfilter2 = Sparse::instance(&shape, &proto2);

        assert!(bloomfilter.merge_inplace(&bloomfilter2).is_ok());
        assert_eq!(*bloomfilter.indices(), [0, 1, 16]);

        assert!(bloomfilter.contains(&bloomfilter2));
    }

    #[test]
    fn filter_merge_inplace_test_sparse_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let mut bloomfilter = Sparse::instance(&shape, &proto);
        let proto2 = SimpleProto::from_u32(0x100);
        let bloomfilter2 = Simple::instance(&shape, &proto2);

        assert!(bloomfilter.merge_inplace(&bloomfilter2).is_ok());
        assert_eq!(*bloomfilter.indices(), [0, 1, 16]);

        assert!(bloomfilter.contains(&bloomfilter2));
    }

    #[test]
    fn filter_merge_inplace_test_sparse_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let mut bloomfilter = Sparse::instance(&shape, &proto);
        let proto2 = SimpleProto::from_u32(0x100);
        let bloomfilter2 = Sparse::instance(&shape, &proto2);

        assert!(bloomfilter.merge_inplace(&bloomfilter2).is_ok());
        assert_eq!(*bloomfilter.indices(), [0, 1, 16]);

        assert!(bloomfilter.contains(&bloomfilter2));
    }

    #[test]
    fn filter_merge_test_simple_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let bloomfilter = Simple::instance(&shape, &proto);
        let proto2 = SimpleProto::from_u32(0x100);
        let bloomfilter2 = Simple::instance(&shape, &proto2);

        let bloomfilter3 = bloomfilter.merge(&bloomfilter2).unwrap();
        assert_eq!(*bloomfilter3.indices(), [0, 1, 16]);
        assert!(bloomfilter3.contains(&bloomfilter));
        assert!(bloomfilter3.contains(&bloomfilter2));
    }

    #[test]
    fn filter_merge_test_simple_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let bloomfilter = Simple::instance(&shape, &proto);
        let proto2 = SimpleProto::from_u32(0x100);
        let bloomfilter2 = Sparse::instance(&shape, &proto2);

        let bloomfilter3 = bloomfilter.merge(&bloomfilter2).unwrap();
        assert_eq!(*bloomfilter3.indices(), [0, 1, 16]);
        assert!(bloomfilter3.contains(&bloomfilter));
        assert!(bloomfilter3.contains(&bloomfilter2));
    }

    #[test]
    fn filter_merge_test_sparse_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let bloomfilter = Sparse::instance(&shape, &proto);
        let proto2 = SimpleProto::from_u32(0x100);
        let bloomfilter2 = Simple::instance(&shape, &proto2);

        let bloomfilter3 = bloomfilter.merge(&bloomfilter2).unwrap();
        assert_eq!(*bloomfilter3.indices(), [0, 1, 16]);
        assert!(bloomfilter3.contains(&bloomfilter));
        assert!(bloomfilter3.contains(&bloomfilter2));
    }

    #[test]
    fn filter_merge_test_sparse_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let bloomfilter = Sparse::instance(&shape, &proto);
        let proto2 = SimpleProto::from_u32(0x100);
        let bloomfilter2 = Sparse::instance(&shape, &proto2);

        let bloomfilter3 = bloomfilter.merge(&bloomfilter2).unwrap();
        assert_eq!(*bloomfilter3.indices(), [0, 1, 16]);
        assert!(bloomfilter3.contains(&bloomfilter));
        assert!(bloomfilter3.contains(&bloomfilter2));
    }

    #[test]
    fn filter_merge_test_simple_by_proto() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let bloomfilter = Simple::instance(&shape, &proto);
        let proto2 = SimpleProto::from_u32(0x100);

        let bloomfilter3 = bloomfilter.merge_proto(&proto2).unwrap();
        assert_eq!(*bloomfilter3.indices(), [0, 1, 16]);
        assert!(bloomfilter3.contains(&bloomfilter));
    }

    #[test]
    fn filter_merge_test_sparse_by_proto() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::from_u32(1);
        let bloomfilter = Sparse::instance(&shape, &proto);
        let proto2 = SimpleProto::from_u32(0x100);

        let bloomfilter3 = bloomfilter.merge_proto(&proto2).unwrap();
        assert_eq!(*bloomfilter3.indices(), [0, 1, 16]);
        assert!(bloomfilter3.contains(&bloomfilter));
    }
}

pub struct BloomFilterFactory {}

impl BloomFilterFactory {
    /// Materializes the proto as a standard Sparse or Simple instance.
    ///
    /// # Arguments
    /// * `shape` - the shape of the filter to materialize
    /// * `proto` - A proto Bloom filter to materialize from.
    fn materialize(shape: &Shape, proto: &dyn Proto) -> BloomFilterType {
        if proto.size() * shape.k < shape.m {
            Sparse::instance(shape, proto)
        } else {
            Simple::instance(shape, proto)
        }
    }
}
