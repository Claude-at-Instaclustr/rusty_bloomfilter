pub mod proto;

use bitvector::*;
use std::vec::Vec;
use std::rc::Rc;
use crate::bloomfilter::proto::Proto;

type BloomFilterType = Box<dyn BloomFilter>;

/// The interface between the internal and external representations of the contents of a Bloom filter.
/// All bloom filters must be able to produce an iterator on bit buckets that map to a bit wise
/// representation of the filter such that for any bit `n`:
///  * `n / 32` yeilds the bit bucket containing the bit.
///  * `1 << (n % 32)` produces the mask to test the bit.
///

/// The traits that all BloomFilters must share
pub trait BloomFilter {

    /// Returns true if it is more efficient to access the bits by index rather than by BitVector.
    /// This is generally the case if hamming_value < shape.m % 8
    fn is_sparse(&self) -> bool {
        self.hamming_value() < self.shape().number_of_buckets()
    }

    fn vector(&self) -> Rc::<BitVector>;

    /// return a list of the bits that are turned on.
    fn indicies(&self) -> Rc::<Vec::<usize>>;

    /// return the shape of the filter
    fn shape(&self) -> &Shape;

    /// Gets the hamming value (the number of bits  turnd on).
    fn hamming_value(&self) -> usize;

    /// Determines if this filter contains the other filter.
    ///
    /// #Errors
    /// Returns an error if the shapes do not match.
    ///
    fn contains(&self, other: &BloomFilterType) -> Result<bool, &str> {
        if self.shape().equivalent_to(other.shape()) {
            print!( "Shapes do not match {:#?} {:#?}", self.shape(), other.shape());
            return Err("Shapes do not match");
        }
        // this is counter intuitive.  We skip all the matches and find the first
        // non matching BitBucket.  If after this there next() return None then the
        // the filters match.
        if self.is_sparse() {
                print!( "self is sparse, ");
                if other.is_sparse()
                {
                    print!( " other is sparse\n");
                    if other.indicies().len() > self.indicies().len()
                    {
                        return Ok(false);
                    }
                    return Ok(other.indicies().iter().all( |s| self.indicies().binary_search( &s ).is_ok()));
                }
                else {
                    print!( "other is not sparse\n");
                    let x : & BitVector = &other.vector();
                    return Ok(x.into_iter().all( |s| self.indicies().binary_search( &s ).is_ok()));
                }
            } else {
                print!( "self is not sparse, ");
                if other.is_sparse()
                {
                    print!( "other is sparse\n");
                    return Ok( other.indicies().iter().all( |s| self.vector().contains( *s ) ));
                } else {
                    print!( "other is not sparse\n");
                    let v1 = self.vector();
                    print!( "v1 = {} {:#?}\n", v1.len(), &v1 );
                    let v2 = other.vector();
                    print!( "v2 = {} {:#?}\n", v2.len(), &v2 );
                    if v1.len() < v2.len()
                    {
                        return Ok( false );
                    }
                    return Ok( v1.iter()
                    .zip( v2.iter())
                    .all( |(s1,s2)| (s1 & s2) == s2 ));
                }

            }
    }

    //
    //  SIZE METHODS
    //

    /// Estimates the number of items in this filter
    fn estimate_n(&self) -> f32 {
        return self.shape().estimate_n(self.hamming_value());
    }

    /// Estimates the number of items in the union of the two filters.
    fn estimate_union(&self, other: &BloomFilterType) -> Result<f32, &str> {
        if self.shape().equivalent_to(other.shape()) {
            return Err("Shapes do not match");
        }
        return Ok( self.shape().estimate_n(self.vector().union( &other.vector() ).len()));
    }

    // Estimates the number of items in the intersection of the two filters.
    fn estimate_intersection(&self, other: &BloomFilterType) -> Result<f32, &str> {
        if self.shape().equivalent_to(other.shape()) {
            return Err("Shapes do not match");
        }
        return Ok(self.estimate_n() + other.estimate_n() - self.estimate_union(other).unwrap());
    }
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
    ///  Gets the
    fn number_of_buckets(&self) -> usize {
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
    shape: Shape,
    buffer: Rc::<BitVector>,
}

impl Simple {

    pub fn empty_instance(shape : &Shape) -> BloomFilterType {
        Box::new(Simple {
            shape: shape.clone(),
            buffer: Rc::new(BitVector::new( shape.number_of_buckets() )),
        })
    }

    pub fn instance(shape: &Shape, proto: & dyn Proto ) -> BloomFilterType {
        let mut v = BitVector::new( shape.number_of_buckets() );
        proto.bits(&shape).iter().for_each( |s| {v.insert( *s );} );
        let simple = Simple {
            shape: shape.clone(),
            buffer: Rc::new(v),
        };
        return Box::new(simple);
    }

}

impl BloomFilter for Simple {

    fn is_sparse(&self) -> bool { false }

    /// return the filter as a BitVector.
    fn vector(&self) -> Rc::<BitVector> {
        self.buffer.clone()
    }

    /// return a list of the bits that are turned on.
    fn indicies(&self) -> Rc::<Vec::<usize>> {
        let x : & BitVector = &self.buffer;
        Rc::new(x.into_iter().collect())
    }

    fn shape(&self) -> &Shape {
        return &self.shape;
    }

    fn hamming_value(&self) -> usize {
        return self.buffer.len();
    }
}

/// A bloom filter that stores the bits in a BitVector
#[derive(Debug)]
pub struct Sparse {
    shape: Shape,
    buffer: Rc::<Vec::<usize>>,
}

impl Sparse {

    pub fn empty_instance(shape : &Shape) -> BloomFilterType {
        Box::new(Sparse {
            shape: shape.clone(),
            buffer: Rc::new(Vec::new()),
        })
    }

    pub fn instance(shape: &Shape, proto: & dyn Proto ) -> BloomFilterType {
        let v = proto.bits(&shape).iter().map( |x| *x ).collect();
        let sparse = Sparse {
            shape: shape.clone(),
            buffer: Rc::new(v),
        };
        return Box::new(sparse);
    }

}

impl BloomFilter for Sparse {

    fn is_sparse(&self) -> bool { true }

    /// return the filter as a BitVector.
    fn vector(&self) -> Rc::<BitVector> {
        let mut v = BitVector::new( self.shape.number_of_buckets() );
        self.buffer.iter().for_each( |s| {v.insert( *s );} );
        Rc::new( v )
    }

    /// return a list of the bits that are turned on.
    fn indicies(&self) -> Rc<Vec::<usize>> {
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

    use crate::bloomfilter::*;
    use crate::bloomfilter::proto::*;

    #[test]
    fn shape_false_positives() {
        let shape = Shape { m: 134_191, k: 23 };
        assert!(shape.false_positives(4000) - (1.0 / 9_994_297.0) < 0.0000001);
    }

    #[test]
    fn shape_number_of_buckets() {
        let shape = Shape { m: 60, k: 2 };
        assert_eq!(shape.number_of_buckets(), 1 );
        let shape = Shape { m: 120, k: 2 };
        assert_eq!(shape.number_of_buckets(), 2 );
    }

    #[test]
    fn empty_filter() {
        let shape = Shape { m: 60, k: 2 };
        let bloomfilter = Simple::empty_instance( &shape );
        assert_eq!( *bloomfilter.indicies(), [] );
        let v = bloomfilter.vector();
        assert_eq!( v.len(), 0 );
        assert_eq!( bloomfilter.hamming_value(), 0 );
        assert!(bloomfilter.estimate_n() < 0.05);
        // filter always contains itself
        assert!(bloomfilter.contains(&bloomfilter).unwrap());
    }

    #[test]
    fn filter_build_correct() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::new(1);
        let bloomfilter = Simple::instance( &shape, &proto );
        assert_eq!( *bloomfilter.indicies(), [0,1] );
        let v = bloomfilter.vector();
        assert_eq!( v.len(), 2 );
        assert!( v.contains( 0 ) );
        assert!( v.contains( 1 ) );
        assert_eq!( bloomfilter.hamming_value(), 2 );
        assert!(bloomfilter.estimate_n() - 1.0 < 0.05);
        // filter always contains itself
        assert!(bloomfilter.contains(&bloomfilter).unwrap());
        let empty_filter: BloomFilterType = Box::new(Simple {
            shape: shape.clone(),
            buffer: Rc::new(BitVector::new( shape.number_of_buckets() )),
        });
        // a filter always contains the empty filter
        assert!(bloomfilter.contains(&empty_filter).unwrap());
        // an empty filter never contains a populated filter
        print!( "{:#?}", empty_filter.contains(&bloomfilter) );
        assert!(!empty_filter.contains(&bloomfilter).unwrap());
        // an empty filter always contains itself
        assert!(empty_filter.contains(&empty_filter).unwrap());
    }

    #[test]
    fn shape_used_multiple_times() {
        let shape = Shape { m: 60, k: 2 };
        let proto = SimpleProto::new(1);
        let bloomfilter = Simple::instance( &shape, &proto );
        let bloomfilter2 = Simple::instance( &shape, &proto );
        assert_eq!( *bloomfilter.indicies(), [0,1] );
        assert_eq!( *bloomfilter2.indicies(), [0,1] );
        assert!( bloomfilter.contains( &bloomfilter2 ).unwrap() );
    }

    #[test]
    fn proto_collection() {
        let shape = Shape { m: 60, k: 2 };
        // this proto will turn on the left 'k' most bits
        let proto = SimpleProto::new(1);
        // this proto will turn on the every other bit for a total of 'k' bits
        let proto2 = SimpleProto::new(2);
        let mut collection = ProtoCollection::new();
        collection.add(Box::new(proto));
        collection.add(Box::new(proto2));
        assert_eq!(collection.count(), 2);
        let bloomfilter = Simple::instance( &shape, &collection );
        assert_eq!( *bloomfilter.indicies(), [0,1,2]);

        //
        // test collection containing a collection
        //

        // this should yeild bits 0 and 16 (65536)
        let proto3 = SimpleProto::new(0x100);
        let mut collection2 = ProtoCollection::new();
        collection2.add(Box::new(collection));
        collection2.add(Box::new(proto3));
        let bloomfilter2 = Simple::instance( &shape, &collection2 );
        assert_eq!( *bloomfilter2.indicies(), [0,1,2,16]);
    }
}
