//pub mod hasher;
pub mod bitmap;
pub mod bitmap_producer;
pub mod hasher;
pub mod index_producer;

use std::collections::HashSet;
use std::vec::Vec;
use crate::bloomfilter::bitmap::BitMap;
use crate::bloomfilter::bitmap_producer::BitMapProducer;
use crate::bloomfilter::hasher::HasherType;
use crate::bloomfilter::index_producer::{IndexProducer, IndexProducerType};


type BloomFilterType = Box<dyn BloomFilter>;

/// The interface between the internal and external representations of the contents of a Bloom filter.
/// All bloom filters must be able to produce an iterator on bit buckets that map to a bit wise
/// representation of the filter such that for any bit `n`:
///  * `n / 32` yeilds the bit bucket containing the bit.
///  * `1 << (n % 32)` produces the mask to test the bit.
///

/// The traits that all BloomFilters must share
pub trait BloomFilter: BitMapProducer + IndexProducer {


    /// Gets the cardinality (number of enabled bits) of this Bloom filter.
    /// This is also known as the Hamming value or Hamming number.</p>
    ///
    fn cardinality(&self) -> usize;

    fn contains_bitmaps( &self, other : &[u64]) -> bool;
    fn contains_hasher(&self, other: &HasherType) -> bool;
    fn contains_indices( &self, other : &HashSet<i32>) -> bool;
    
       // Estimates the number of items in the intersection of the two filters.
    fn estimate_intersection(&self, other: &BloomFilterType) -> f32 {
        return self.estimate_n() + other.estimate_n() - self.estimate_union(other);
    }
    
     /// Estimates the number of items in this filter
    fn estimate_n(&self) -> f32 {
        return self.shape().estimate_n(self.cardinality() as usize);
    }
    
    /// Estimates the number of items in the union of the two filters.
    fn estimate_union(&self, other: &BloomFilterType) -> f32 {
        let estimate = match self.is_sparse() && other.is_sparse() {
            true => self.merge_indices( &other.get_indices() ),
            false => self.merge_bitmaps( other.get_bitmaps().as_slice() ),
        };
        estimate.estimate_n()
    }

    /// Returns the shape of the filter
    fn shape(&self) -> &Shape;
    
    fn is_full(&self) -> bool {
    	self.cardinality() == self.shape().m
    }
    
    /// Returns true if it is more efficient to access the bits by index rather than by BitVector.
    /// This is generally the case if hamming_value < shape.m % 8
    fn is_sparse(&self) -> bool;

    /// Aggregates two Bloom filtrs.  If the merge will probably produce a sparse filter the sparse
    /// implementation is used.  Otherwise the Simple implementaiton is used.
    /// If the bloom filter arguments are not of Simple or Sparse they will be treated as though
    /// they were.  (e.g. no special handling for other types)
    fn merge_indices(&self, other: &HashSet<i32>) -> BloomFilterType;

    fn merge_bitmaps(&self, other: &[u64]) -> BloomFilterType;

    fn merge_hasher(&self, hasher: &HasherType) -> BloomFilterType;
    
    // Updates this filter by merging the values from the other filter
    fn merge_bitmaps_in_place(&mut self, other: &[u64]) -> bool;

    fn merge_indices_in_place(&mut self, other: &HashSet<i32>) -> bool;
    fn merge_hasher_in_place(&mut self, hasher: &HasherType) -> bool;
}

/// The shape of the bloom filter.
#[derive(Copy, Clone, Debug)]
pub struct Shape {
    /// number of bits in the filter
    m: usize,
    /// number of functions
    k: usize,
}

impl Shape {

    pub fn is_sparse( &self, bits : usize ) -> bool {
        // if there are 2 bits per bucket then the sparse and the non spare have the
        // same storage requirements.
        self.number_of_buckets() <= (2 * bits) as usize
    }

    ///  Gets the
    pub fn number_of_buckets(&self) -> usize {
        let mut len = self.m / 64;
        if self.m % 64 > 0 {
            len += 1;
        }
        len
    }

    pub fn equivalent_to(&self, other: &Shape) -> bool {
        self.m != other.m || self.k != other.k
    }
    
    /// Calculaes the probability of false positives if `count` objects are added to a filter using this shape.
    pub fn false_positives(&self, count: usize) -> f64 {
        let k = self.k as f64;
        let n = count as f64;
        let m = self.m as f64;
        (1.0 - ((-k * n) / m).exp()).powf(k)
    }

    /// Calculates the estimated number of items in a filter
    /// of this Shape.
    ///
    /// count is the number of bits turned on in the filter.
    pub fn estimate_n(&self, count: usize) -> f32 {
        let c = count as f64;
        let m = self.m as f64;
        let k = self.k as f64;
        let result = -(m / k) * (1.0 - (c / m)).ln();
        result as f32
    }
}

/// A bloom filter that stores the bits in a BitVector
#[derive(Debug)]
pub struct Simple {
    cardinality : usize,
    shape: Shape,
    buffer: Vec<u64>,
}

impl Simple {
    pub fn empty_instance(shape: &Shape) -> BloomFilterType {
        Box::new(Simple {
            cardinality : 0,
            shape: shape.clone(),
            buffer: Vec::with_capacity( shape.number_of_buckets() ),
        })
    }

    pub fn from_hasher(shape: &Shape, hasher: &HasherType) -> BloomFilterType {
        let mut result = Simple::empty_instance(&shape);
        result.merge_hasher_in_place( hasher );
        result
    }

    fn calc_cardinality(&self) -> usize {
        let mut result :usize = 0;
        for word in &self.buffer {
            result += word.count_ones() as usize;
        }
        result
    }

}

impl IndexProducer for Simple {

    fn get_indices(&self) -> HashSet<i32> {
        let mut result : HashSet<i32> = HashSet::with_capacity( self.cardinality );
        for block in  0..self.buffer.len() {

            if self.buffer[block] != 0 {
                let mut mask : u64 = 1;
                for offset in 0..64  {
                    if (self.buffer[block] & mask) == mask {
                        result.insert(((block * 64) + offset) as i32);
                    }
                    mask <<= 1;
                }
            }
        }
        result
    }
}

impl BitMapProducer for Simple {
    fn get_bitmaps(&self) -> Vec<u64> {
        self.buffer.clone()
    }
}

impl BloomFilter for Simple {
    fn cardinality(&self) -> usize {
        self.cardinality
    }

    fn contains_bitmaps( &self, other : &[u64]) -> bool {
        if other.len() > self.buffer.len() {
            return false;
        }

        for idx in 0..other.len() {
            if (self.buffer[idx] & other[idx]) != other[idx] {
                return false;
            }
        }
        true
    }
    fn contains_hasher(&self, other: &HasherType) -> bool {
        let producers : Vec<IndexProducerType> = other.indices( self.shape() );
        for producer in producers {
            if ! self.contains_indices( &producer.get_indices() ) {
                return false;
            }
        }
        return true;
    }

    fn contains_indices( &self, other : &HashSet<i32>) -> bool {
        for idx in other {
            if ! BitMap::contains(&self.buffer, idx ) {
                return false;
            }
        }
        true
    }

    fn shape(&self) -> &Shape {
        return &self.shape;
    }

    fn is_sparse(&self) -> bool {
        false
    }

    fn merge_indices(&self, other : &HashSet<i32> ) -> BloomFilterType {
        let mut simple = Simple {
            cardinality : self.cardinality,
            shape : self.shape.clone(),
            buffer : self.buffer.clone()
        };
        simple.merge_indices_in_place(other);
        Box::new( simple )
    }

    fn merge_bitmaps(&self, other : &[u64] ) -> BloomFilterType {
        let mut simple = Simple {
            cardinality : self.cardinality,
            shape : self.shape.clone(),
            buffer : self.buffer.clone()
        };
        simple.merge_bitmaps_in_place( other );
        Box::new(simple )
    }

    fn merge_hasher(&self, hasher: &HasherType) -> BloomFilterType {
            let mut result = Simple {
                cardinality : self.cardinality,
                shape : self.shape.clone(),
                buffer : self.buffer.clone(),
            };
            result.merge_hasher_in_place( hasher );
            Box::new( result )
    }

    fn merge_bitmaps_in_place(&mut self, other: &[u64]) -> bool {
        for idx in 0..other.len() {
            while self.buffer.len() <= idx {
                self.buffer.push( 0 );
            }
            self.buffer[idx] |= other[idx];
        }
        self.cardinality = self.calc_cardinality();
        true
    }

    fn merge_indices_in_place(&mut self, other: &HashSet<i32>) -> bool {
        for idx in other {
            BitMap::set( &mut self.buffer, idx );
        }
        self.cardinality = self.calc_cardinality();
        true
    }

    fn merge_hasher_in_place(&mut self, other: &HasherType) -> bool {
        let producers : Vec<IndexProducerType> = other.indices( self.shape() );
        for producer in producers {
            if ! self.merge_indices_in_place( &producer.get_indices() ) {
                return false;
            }
        }
        true
    }
}


/// A bloom filter that stores the bits in a BitVector
#[derive(Debug)]
pub struct Sparse {
    shape: Shape,
    buffer: HashSet<i32>,
}

impl Sparse {
    pub fn empty_instance(shape: &Shape) -> BloomFilterType {
        Box::new(Sparse {
            shape: shape.clone(),
            buffer: HashSet::new(),
        })
    }

    pub fn from_hasher(shape: &Shape, hasher: &HasherType) -> BloomFilterType {
        let mut result = Sparse::empty_instance(shape);
        for producer in hasher.indices(shape) {
            result.merge_indices_in_place(&producer.get_indices());
        }
        result
    }

}

impl IndexProducer for Sparse {
    fn get_indices(&self) -> HashSet<i32> {
        self.buffer.clone()
    }
}

impl BitMapProducer for Sparse {
    fn get_bitmaps(&self) -> Vec<u64> {
        let mut bitmaps: Vec<u64> = Vec::with_capacity( self.shape.number_of_buckets());
        for i in self.buffer.iter() {
            BitMap::set(&mut bitmaps, i );
        }
        return bitmaps;
    }
}

impl BloomFilter for Sparse {
    fn cardinality(&self)-> usize {
        self.buffer.len()
    }

    fn contains_hasher(&self, other: &HasherType) -> bool {
        for producer in other.indices( self.shape()) {
            if ! self.contains_indices( &producer.get_indices() ) {
                return false;
            }
        }
        true
    }

    fn contains_bitmaps( &self, other : &[u64]) -> bool {
        for block in  0..other.len() {
            if other[block] != 0 {
                let mut mask : u64 = 1;
                for offset in 0..64  {
                    if (other[block] & mask) == mask {
                        let idx = ((block * 64) + offset) as i32;
                        if !self.buffer.contains(&idx) {
                            return false;
                        }
                    }
                    mask <<= 1;
                }
            }
        }
     true
    }

    fn contains_indices( &self, other : &HashSet<i32>) -> bool {
        for idx in other {
            if ! self.buffer.contains( idx ) {
                return false;
            }
        }
        return true;
    }

    fn shape(&self) -> &Shape {
        return &self.shape;
    }

    fn is_sparse(&self) -> bool {
        true
    }

    fn merge_bitmaps(&self, other : &[u64] ) -> BloomFilterType {
        let mut simple = Simple::empty_instance(self.shape());
        simple.merge_bitmaps_in_place(other);
        simple.merge_indices_in_place( &self.get_indices());
        simple
    }

    fn merge_indices(&self, other : &HashSet<i32> ) -> BloomFilterType {
        let mut simple = Simple::empty_instance(self.shape());
        simple.merge_indices_in_place(other);
        simple.merge_indices_in_place( &self.get_indices());
        simple
    }

    fn merge_indices_in_place(&mut self, other: &HashSet<i32>) -> bool {
        for idx in other {
            self.buffer.insert( *idx );
        }
        true
    }

    fn merge_bitmaps_in_place(&mut self, other: &[u64]) -> bool {
        for block in  0..other.len() {
            if other[block] != 0 {
                let mut mask : u64 = 1;
                for offset in 0..64  {
                    if (other[block] & mask) == mask {
                        if !self.buffer.insert(((block * 64) + offset) as i32) {
                            return false;
                        }
                    }
                    mask <<= 1;
                }
            }
        }
        true
    }

    fn merge_hasher(&self, hasher: &HasherType) -> BloomFilterType  {
        let mut result : BloomFilterType = if self.shape().is_sparse(
            ((hasher.size() as usize)*self.shape().k) + self.cardinality() ){
            // create a sparse one
            print!("Creating sparse BloomFilter\n");
            Sparse::empty_instance( self.shape())
        } else {
            // create a simple one
            print!("Creating simple BloomFilter\n");
            Simple::empty_instance( self.shape())
        };
        result.merge_indices_in_place( &self.get_indices() );
        result.merge_hasher_in_place( hasher );
        result
    }

    fn merge_hasher_in_place(&mut self, other: &HasherType) -> bool {
        let producers : Vec<IndexProducerType> = other.indices( self.shape() );
        for producer in producers {
            if ! self.merge_indices_in_place( &producer.get_indices() ) {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use crate::bloomfilter::*;
    use crate::bloomfilter::hasher::SimpleHasher;

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
        assert_eq!(bloomfilter.get_indices().len(), 0);
        assert_eq!(bloomfilter.get_bitmaps().len(), 0);
        assert_eq!(bloomfilter.cardinality(), 0);
        assert!(bloomfilter.estimate_n() < 0.05);
        // filter always contains itself
        assert!(bloomfilter.contains_indices(&bloomfilter.get_indices()));
        assert!(bloomfilter.contains_bitmaps(&bloomfilter.get_bitmaps()));
    }

    #[test]
    fn filter_build_correct() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let bloomfilter = Simple::from_hasher(&shape, &hasher);
        let v = bloomfilter.get_indices();
        assert_eq!(v.len(), 2);
        assert!(v.contains(&0));
        assert!(v.contains(&1));
        assert_eq!(bloomfilter.cardinality(), 2);
        assert!(bloomfilter.estimate_n() - 1.0 < 0.05);
        // filter always contains itself
        assert!(bloomfilter.contains_hasher(&hasher));
        assert!(bloomfilter.contains_indices(&bloomfilter.get_indices()));
        assert!(bloomfilter.contains_bitmaps(&bloomfilter.get_bitmaps()));
        let empty_filter: BloomFilterType = Simple::empty_instance(&shape );
        // a filter always contains the empty filter
        assert!(bloomfilter.contains_indices(&empty_filter.get_indices()));
        assert!(bloomfilter.contains_bitmaps(&empty_filter.get_bitmaps()));
        // an empty filter never contains a populated filter
        assert!(! empty_filter.contains_indices(&bloomfilter.get_indices()));
        assert!(! empty_filter.contains_bitmaps(&bloomfilter.get_bitmaps()));
        assert!(! empty_filter.contains_hasher(&hasher));
        // an empty filter always contains itself
        assert!(empty_filter.contains_indices(&empty_filter.get_indices()));
        assert!(empty_filter.contains_bitmaps(&empty_filter.get_bitmaps()));
    }

    fn contains_test(bloomfilter : BloomFilterType, bloomfilter2 : BloomFilterType) {

        let simple_empty_filter: BloomFilterType = Simple::empty_instance( bloomfilter.shape() );
        let sparse_empty_filter: BloomFilterType = Sparse::empty_instance( bloomfilter.shape() );

        // filter always contains itself
        assert!(bloomfilter.contains_indices(&bloomfilter.get_indices()));
        assert!(bloomfilter.contains_bitmaps(&bloomfilter.get_bitmaps()));
        // a filter always contains the empty filter
        assert!(bloomfilter.contains_indices(&simple_empty_filter.get_indices()));
        assert!(bloomfilter.contains_bitmaps(&simple_empty_filter.get_bitmaps()));
        assert!(bloomfilter.contains_indices(&sparse_empty_filter.get_indices()));
        assert!(bloomfilter.contains_bitmaps(&sparse_empty_filter.get_bitmaps()));
        // an empty filter never contains a populated filter
        assert!(!simple_empty_filter.contains_indices(&bloomfilter.get_indices()));
        assert!(!simple_empty_filter.contains_bitmaps(&bloomfilter.get_bitmaps()));
        assert!(!sparse_empty_filter.contains_indices(&bloomfilter.get_indices()));
        assert!(!sparse_empty_filter.contains_bitmaps(&bloomfilter.get_bitmaps()));
        // an empty filter always contains itself
        assert!(simple_empty_filter.contains_indices(&simple_empty_filter.get_indices()));
        assert!(simple_empty_filter.contains_bitmaps(&simple_empty_filter.get_bitmaps()));
        assert!(simple_empty_filter.contains_indices(&sparse_empty_filter.get_indices()));
        assert!(simple_empty_filter.contains_bitmaps(&sparse_empty_filter.get_bitmaps()));
        assert!(sparse_empty_filter.contains_indices(&simple_empty_filter.get_indices()));
        assert!(sparse_empty_filter.contains_bitmaps(&simple_empty_filter.get_bitmaps()));
        assert!(sparse_empty_filter.contains_indices(&sparse_empty_filter.get_indices()));
        assert!(sparse_empty_filter.contains_bitmaps(&sparse_empty_filter.get_bitmaps()));

        // bloom filter2 contains bloom filter 1
        assert!(bloomfilter2.contains_indices(&bloomfilter.get_indices()));
        assert!(bloomfilter2.contains_bitmaps(&bloomfilter.get_bitmaps()));
        // bloom filter1 does not contain bloom filter 2
        assert!(!bloomfilter.contains_indices(&bloomfilter2.get_indices()));
        assert!(!bloomfilter.contains_bitmaps(&bloomfilter2.get_bitmaps()));
    }

    #[test]
    fn contains_test_simple_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0, 1);
        let hasher2 = SimpleHasher::new(0, 5);
        let bloomfilter = Simple::from_hasher(&shape, &hasher);

        let mut bloomfilter2 = Simple::from_hasher(&shape, &hasher);
        assert!(bloomfilter2.merge_hasher_in_place(&hasher2));

        contains_test( bloomfilter, bloomfilter2);
    }

    #[test]
    fn contains_test_simple_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let hasher2 = SimpleHasher::new(0,5);
        let bloomfilter = Simple::from_hasher(&shape, &hasher);

        let mut bloomfilter2 = Sparse::from_hasher(&shape, &hasher);
        assert!(bloomfilter2.merge_hasher_in_place( &hasher2 ));

        contains_test( bloomfilter, bloomfilter2);
    }

    #[test]
    fn contains_test_sparse_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let hasher2 = SimpleHasher::new(0,5);
        let bloomfilter = Sparse::from_hasher(&shape, &hasher);

        let mut bloomfilter2 = Simple::from_hasher(&shape, &hasher);
        assert!(bloomfilter2.merge_hasher_in_place(&hasher2));

        contains_test( bloomfilter, bloomfilter2);

    }

    #[test]
    fn contains_test_sparse_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let hasher2 = SimpleHasher::new(0,5);
        let bloomfilter = Sparse::from_hasher(&shape, &hasher);

        let mut bloomfilter2 = Sparse::from_hasher(&shape, &hasher);
        assert!(bloomfilter2.merge_hasher_in_place( &hasher2 ));

        contains_test( bloomfilter, bloomfilter2);

    }

    #[test]
    fn shape_used_multiple_times() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let bloomfilter = Simple::from_hasher(&shape, &hasher);
        let bloomfilter2 = Simple::from_hasher(&shape, &hasher);
        assert_eq!(bloomfilter.get_indices().len(), 2);
        assert!( bloomfilter.get_indices().contains( &0 ));
        assert!( bloomfilter.get_indices().contains( &1 ));

        assert_eq!(bloomfilter2.get_indices().len(), 2);
        assert!( bloomfilter2.get_indices().contains( &0 ));
        assert!( bloomfilter2.get_indices().contains( &1 ));

        assert!(bloomfilter.contains_bitmaps(&bloomfilter2.get_bitmaps()));
    }

    #[test]
    fn filter_merge_inplace_test_simple_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let mut bloomfilter = Simple::from_hasher(&shape, &hasher);
        let hasher2 = SimpleHasher::new(0,0x100);
        let bloomfilter2 = Simple::from_hasher(&shape, &hasher2);

        assert!(bloomfilter.merge_hasher_in_place(&hasher2));
        assert_eq!(bloomfilter.get_indices().len(), 3);
        assert!( bloomfilter.get_indices().contains( &0 ));
        assert!( bloomfilter.get_indices().contains( &1 ));
        assert!( bloomfilter.get_indices().contains( &16 ));

        assert!(bloomfilter.contains_bitmaps(&bloomfilter2.get_bitmaps()));
    }

    #[test]
    fn filter_merge_inplace_test_simple_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let mut bloomfilter = Simple::from_hasher(&shape, &hasher);
        let hasher2 = SimpleHasher::new(0, 0x100);
        let bloomfilter2 = Sparse::from_hasher(&shape, &hasher2);

        assert!(bloomfilter.merge_hasher_in_place(&hasher2));
        assert_eq!(bloomfilter.get_indices().len(), 3);
        assert!( bloomfilter.get_indices().contains( &0 ));
        assert!( bloomfilter.get_indices().contains( &1 ));
        assert!( bloomfilter.get_indices().contains( &16 ));

        assert!(bloomfilter.contains_bitmaps(&bloomfilter2.get_bitmaps()));
    }

    #[test]
    fn filter_merge_inplace_test_sparse_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let mut bloomfilter = Sparse::from_hasher(&shape, &hasher);
        let hasher2 = SimpleHasher::new(0, 0x100);
        let bloomfilter2 = Simple::from_hasher(&shape, &hasher2);

        assert!(bloomfilter.merge_hasher_in_place(&hasher2));
        assert_eq!(bloomfilter.get_indices().len(), 3);
        assert!( bloomfilter.get_indices().contains( &0 ));
        assert!( bloomfilter.get_indices().contains( &1 ));
        assert!( bloomfilter.get_indices().contains( &16 ));

        assert!(bloomfilter.contains_bitmaps(&bloomfilter2.get_bitmaps()));
    }

    #[test]
    fn filter_merge_inplace_test_sparse_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let mut bloomfilter = Sparse::from_hasher(&shape, &hasher);
        let hasher2 = SimpleHasher::new(0, 0x100);
        let bloomfilter2 = Sparse::from_hasher(&shape, &hasher2);

        assert!(bloomfilter.merge_hasher_in_place(&hasher2));
        assert_eq!(bloomfilter.get_indices().len(), 3);
        assert!( bloomfilter.get_indices().contains( &0 ));
        assert!( bloomfilter.get_indices().contains( &1 ));
        assert!( bloomfilter.get_indices().contains( &16 ));

        assert!(bloomfilter.contains_bitmaps(&bloomfilter2.get_bitmaps()));
    }

    #[test]
    fn filter_merge_test_simple_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let bloomfilter = Simple::from_hasher(&shape, &hasher);
        let hasher2 = SimpleHasher::new(0, 0x100);
        let bloomfilter2 = Simple::from_hasher(&shape, &hasher2);

        let bloomfilter3 = bloomfilter.merge_bitmaps(&bloomfilter2.get_bitmaps());
        assert_eq!(bloomfilter3.get_indices().len(), 3);
        assert!( bloomfilter3.get_indices().contains( &0 ));
        assert!( bloomfilter3.get_indices().contains( &1 ));
        assert!( bloomfilter3.get_indices().contains( &16 ));
        assert!(bloomfilter3.contains_bitmaps(&bloomfilter.get_bitmaps()));
        assert!(bloomfilter3.contains_bitmaps(&bloomfilter2.get_bitmaps()));
    }

    #[test]
    fn filter_merge_test_simple_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let bloomfilter = Simple::from_hasher(&shape, &hasher);
        let hasher2 = SimpleHasher::new(0, 0x100);
        let bloomfilter2 = Sparse::from_hasher(&shape, &hasher2);

        let bloomfilter3 = bloomfilter.merge_indices(&bloomfilter2.get_indices());
        assert_eq!(bloomfilter3.get_indices().len(), 3);
        assert!( bloomfilter3.get_indices().contains( &0 ));
        assert!( bloomfilter3.get_indices().contains( &1 ));
        assert!( bloomfilter3.get_indices().contains( &16 ));
        assert!(bloomfilter3.contains_bitmaps(&bloomfilter.get_bitmaps()));
        assert!(bloomfilter3.contains_bitmaps(&bloomfilter2.get_bitmaps()));
    }

    #[test]
    fn filter_merge_test_sparse_by_simple() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let bloomfilter = Sparse::from_hasher(&shape, &hasher);
        let hasher2 = SimpleHasher::new(0, 0x100);
        let bloomfilter2 = Simple::from_hasher(&shape, &hasher2);

        let bloomfilter3 = bloomfilter.merge_indices(&bloomfilter2.get_indices());
        let bloomfilter4 = bloomfilter.merge_bitmaps(&bloomfilter2.get_bitmaps());
        assert_eq!(bloomfilter3.get_indices().len(), 3);
        assert!( bloomfilter3.get_indices().contains( &0 ));
        assert!( bloomfilter3.get_indices().contains( &1 ));
        assert!( bloomfilter3.get_indices().contains( &16 ));
        assert!(bloomfilter3.contains_bitmaps(&bloomfilter.get_bitmaps()));
        assert!(bloomfilter3.contains_bitmaps(&bloomfilter2.get_bitmaps()));
        assert!(bloomfilter4.contains_bitmaps(&bloomfilter.get_bitmaps()));
        assert!(bloomfilter4.contains_bitmaps(&bloomfilter2.get_bitmaps()));
    }

    #[test]
    fn filter_merge_test_sparse_by_sparse() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let bloomfilter = Sparse::from_hasher(&shape, &hasher);
        let hasher2 = SimpleHasher::new(0, 0x100);
        let bloomfilter2 = Sparse::from_hasher(&shape, &hasher2);

        let bloomfilter3 = bloomfilter.merge_indices(&bloomfilter2.get_indices());
        let bloomfilter4 = bloomfilter.merge_bitmaps(&bloomfilter2.get_bitmaps());
        assert_eq!(bloomfilter3.get_indices().len(), 3);
        assert!( bloomfilter3.get_indices().contains( &0 ));
        assert!( bloomfilter3.get_indices().contains( &1 ));
        assert!( bloomfilter3.get_indices().contains( &16 ));
        assert!(bloomfilter3.contains_bitmaps(&bloomfilter.get_bitmaps()));
        assert!(bloomfilter3.contains_bitmaps(&bloomfilter2.get_bitmaps()));
        assert!(bloomfilter4.contains_bitmaps(&bloomfilter.get_bitmaps()));
        assert!(bloomfilter4.contains_bitmaps(&bloomfilter2.get_bitmaps()));
    }

    #[test]
    fn filter_merge_test_simple_by_hasher() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let bloomfilter = Simple::from_hasher(&shape, &hasher);
        let hasher2 = SimpleHasher::new(0, 0x100);

        let bloomfilter3 = bloomfilter.merge_hasher(&hasher2);
        assert_eq!(bloomfilter3.get_indices().len(), 3);
        assert!( bloomfilter3.get_indices().contains( &0 ));
        assert!( bloomfilter3.get_indices().contains( &1 ));
        assert!( bloomfilter3.get_indices().contains( &16 ));

        assert!(bloomfilter3.contains_bitmaps(&bloomfilter.get_bitmaps()));
        assert!(bloomfilter3.contains_hasher(&hasher));
        assert!(bloomfilter3.contains_hasher(&hasher2));
    }

    #[test]
    fn filter_merge_test_sparse_by_hasher() {
        let shape = Shape { m: 60, k: 2 };
        let hasher = SimpleHasher::new(0,1);
        let bloomfilter = Sparse::from_hasher(&shape, &hasher);
        let hasher2 = SimpleHasher::new(0, 0x100);

        let bloomfilter3 = bloomfilter.merge_hasher(&hasher2);
        assert_eq!(bloomfilter3.get_indices().len(), 3);
        assert!( bloomfilter3.get_indices().contains( &0 ));
        assert!( bloomfilter3.get_indices().contains( &1 ));
        assert!( bloomfilter3.get_indices().contains( &16 ));

        assert!(bloomfilter3.contains_bitmaps(&bloomfilter.get_bitmaps()));
        assert!(bloomfilter3.contains_hasher(&hasher));
        assert!(bloomfilter3.contains_hasher(&hasher2));
    }
}
