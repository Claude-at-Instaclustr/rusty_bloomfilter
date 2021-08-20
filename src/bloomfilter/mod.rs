use std::convert::TryInto;
use std::vec::Vec;

type BloomFilterType = Box<dyn BloomFilter>;

/// The traits that all BloomFilters must share
pub trait BloomFilter {
	/// Turn on the specified bit.
	///
	/// #Panics
	/// Will panic if the selected bit is >=n `shape.m`
	///
	fn on( &mut self, bit : u32 );

	/// Turn off the specified bit.
	///
	/// #Panics
	/// Will panic if the selected bit is >=n `shape.m`
	///
	fn off( &mut self, bit : u32 );

	/// return the filter as a vector of u32 bit vectors.
	fn as_vec( &self ) -> Box<&Vec::<u32>>;

	/// return the shape of the filter
	fn shape( &self ) -> &Shape;

	/// Determines if this filter contains the other filter.
	///
	/// #Errors
	/// Returns an error if the shapes do not match.
	///
	fn contains( &self, other : &BloomFilterType ) -> Result<bool, &str>
	{
		if self.shape().equivalent_to( other.shape() ) {
			return Err( "Shapes do not match" )
		}
		let self_buff = self.as_vec();
		let other_buff = self.as_vec();
		for i in 0..self_buff.len() {
			if self_buff[i] & other_buff[i] != other_buff[i] {
				return Ok(false);
			}
		}
		Ok(true)
	}

	/// Gets the hamming value (the number of bits  turnd on).
	fn hamming_value( &self ) -> u32;


	/// Estimates the number of items in this filter
	fn estimate_n( &self ) -> f32 {
		return self.shape().estimate_n( self.hamming_value() );
	}

	/// Estimates the number of items in the union of the two filters.
	fn estimate_union( &self,  other : &BloomFilterType ) -> Result<f32, &str> {
		if self.shape().equivalent_to( other.shape() ) {
			return Err( "Shapes do not match" )
		}
		let self_buff = self.as_vec();
		let other_buff = self.as_vec();
		let mut count = 0;
		for i in 0..self_buff.len() {
			let x = self_buff[i] | other_buff[i];
			count += x.count_ones();
		}
		return Ok( self.shape().estimate_n( count ) );
	}

	// estimate the number of items in the intersection of the two filters.
	fn estimate_intersection( &self, other : &BloomFilterType ) -> Result<f32, &str> {
		if self.shape().equivalent_to( other.shape() ) {
			return Err( "Shapes do not match" )
		}
		return Ok(self.estimate_n() + other.estimate_n() - self.estimate_union( other ).unwrap());
	}

}

/// A prototype Bloom filer.
///
/// Proto bloom filters can build filters of any shape.
/// Proto filters process the hash function to build the filter based on shape.
pub trait Proto {
	/// Create a bloom filter from the shape.
	/// The specific type of the filter is implementation specific.
	fn build( &self, shape : &Shape ) -> BloomFilterType;

	/// Add this proto to a filter.
	fn add_to( &self,  filter : &mut BloomFilterType);
}

/// A Proto implementation that is a collection of Protos.
///
/// This represents multiple objects to be placed in the filter.
pub struct ProtoCollection {
	inner: Vec<Box<dyn Proto>>,
}

/// A Proto implementation that is a single object to be placed in the filter.
pub struct SimpleProto {
	start : u32,
	incr : u32,
}

/// The shape of the bloom filter.
#[derive(Copy, Clone)]
pub struct Shape {
	/// number of bits in the filter
	m : u32,
	/// number of functions
	k : u32,
}


/// The simple prototype
impl SimpleProto {
	pub fn new( full : u64 ) -> SimpleProto {
		let start =   (full >> 32) as u32 ;
		let incr = full  as u32;

		return SimpleProto{ start, incr };
	}
}

impl Proto for SimpleProto {
	/// Build a Simple Bloomfilter from the prototype
	fn build( &self, shape : &Shape ) -> BloomFilterType {
		let mut filter =  Simple::instance( shape );
		self.add_to( &mut filter );
		return filter;
	}

	fn add_to( &self,  filter : &mut BloomFilterType ) {
		let mut accumulator = self.start;
		for _i in 0..filter.shape().k {
			filter.on( accumulator % filter.shape().m );
			accumulator = accumulator + self.incr;
		}
	}
}

impl ProtoCollection {
	/// Creates an empty collection
	pub fn new() -> ProtoCollection {
		return ProtoCollection{ inner : Vec::new() }
	}

	/// Add prototypes to the collection
	pub fn add(&mut self, proto : Box<dyn Proto> ) {
		self.inner.push( proto );
	}

	/// Gets the number of Protos in the collection
	pub fn len(&self) -> usize {
		return self.inner.len();
	}
}

impl Proto for ProtoCollection {
	/// create a Simple filter with the specified shape from this proto.
	fn build( &self, shape : &Shape ) -> BloomFilterType {
		let mut filter  =  Simple::instance( shape );
		self.add_to( &mut filter );
		return filter
	}

	fn add_to( &self, filter : &mut BloomFilterType) {
		for v in &self.inner {
			v.add_to(  filter )
		}
	}
}

impl Shape {

	pub fn equivalent_to( &self, other : &Shape ) -> bool {
		return self.m != other.m || self.k != other.k;
	}
	/// Calculaes the probability of false positives if `count` objects are added to a filter using this shape.
	pub fn false_positives( &self, count : u32 ) -> f64 {
		let k = self.k as f64;
		let n = count as f64;
		let m = self.m as f64;
		return ( 1.0-((-k*n)/m).exp()).powf(k);
	}

	/// Calculates the estimated number of items in a filter
	/// of this Shape.
	///
	/// count is the number of bits turned on in the filter.
	pub fn estimate_n( &self, count : u32 ) -> f32 {
		let c = count as f64;
		let m = self.m as f64;
		let k = self.k as f64;
		let result = -(m/k) * (1.0 - (c/m)).ln();
		return result as f32;
	}
}



/// A bloom filter that stores the bits in an array of u32 bit buffers
pub struct Simple {
	shape : Shape,
  	buffer : Vec::<u32>,
}

impl Simple {

	pub fn instance( shape : &Shape ) -> BloomFilterType {
		return Box::new( Simple::new( shape ));
	}

	/// Construct an empty Simple Bloom filter
	pub fn new( shape : &Shape ) -> Simple {
		let size = if shape.m % 32 > 0  {1+(shape.m / 32)} else {shape.m / 32};
		let buff = vec![0 as u32 ; size as usize ];
		Simple{ shape: shape.clone(), buffer : buff }
	}

}

impl BloomFilter for Simple {
 	fn as_vec(&self) -> Box<&Vec::<u32>> {
		return Box::new( &self.buffer );
	}

	fn shape(&self) -> &Shape {
		return &self.shape;
	}

	fn on( &mut self, bit : u32 ) {
		if bit >= self.shape.m {
			panic!( "Bit position too large")
		}
		let offset = bit % 32;
		let position : usize = (bit / 32).try_into().unwrap();
		self.buffer[position] |= 1<<offset;
	}

	fn off( &mut self, bit : u32 ) {
		if bit >= self.shape.m {
			panic!( "Bit position too large")
		}
		let offset = bit % 32;
		let position : usize = (bit / 32).try_into().unwrap() ;
		self.buffer[position] &= !(1<<offset);
	}

	fn hamming_value( &self ) -> u32 {
		return self.buffer.iter().map(|i| i.count_ones() ).sum();
	}


}
//impl BloomFilter {
//	pub fn new( shape : &Shape ) -> BloomFilter {
//
//		let size = if shape.m % 32 > 0  {1+(shape.m / 32)} else {shape.m / 32};
//		let buff = vec![0 as u32 ; size as usize ];
//		BloomFilter{ shape: shape.clone(), buffer : buff }
//	}
//
//
//	pub fn on( &mut self, bit : u32 ) {
//		if bit >= self.shape.m {
//			panic!( "Bit position too large")
//		}
//		let offset = bit % 32;
//		let position : usize = (bit / 32).try_into().unwrap();
//		self.buffer[position] |= 1<<offset;
//	}
//
//	pub fn off( &mut self, bit : u32 ) {
//		if bit >= self.shape.m {
//			panic!( "Bit position too large")
//		}
//		let offset = bit % 32;
//		let position : usize = (bit / 32).try_into().unwrap() ;
//		self.buffer[position] &= !(1<<offset);
//	}
//
//	pub fn contains( &self, filter : &BloomFilter ) -> Result<bool, &str>
//	{
//		if self.shape.m != filter.shape.m || self.shape.k != filter.shape.k {
//			return Err( "Shapes do not match" )
//		}
//		for i in 0..self.buffer.len() {
//			if self.buffer[i] & filter.buffer[i] != filter.buffer[i] {
//				return Ok(false);
//			}
//		}
//		Ok(true)
//	}
//
//	pub fn hamming_value( &self ) -> u32 {
//		return self.buffer.iter().map(|i| i.count_ones() ).sum();
//	}
//
//	pub fn estimate_n( &self ) -> f32 {
//		return self.shape().estimate( self.hamming_value() );
//	}
//
//	pub fn estimate_union( &self,  other : &BloomFilter ) -> f32 {
//		let mut count = 0;
//		for i in 0..self.buffer.len() {
//			let x = self.buffer[i] | other.buffer[i];
//			count += x.count_ones();
//		}
//		return self.estimate( count );
//	}
//
//
//	pub fn estimate_intersection( &self, other : &BloomFilter ) -> f32 {
//		return self.estimate_n() + other.estimate_n() - self.estimate_union( other );
//	}
//}

#[cfg(test)]
mod tests {

	use crate::bloomfilter::BloomFilterType;
	use crate::bloomfilter::Proto;
	use crate::bloomfilter::SimpleProto;
	use crate::bloomfilter::ProtoCollection;
	use crate::bloomfilter::Simple;
	use crate::bloomfilter::Shape;

    #[test]
    fn shape_false_positives() {
		let shape = Shape{ m : 134_191, k : 23};
        assert!(shape.false_positives( 4000 )-(1.0/9_994_297.0) < 0.0000001 );
    }

	#[test]
	fn simple_proto_correct() {
		let proto  = SimpleProto::new( 1 );
		assert_eq!( proto.start, 0 );
		assert_eq!( proto.incr, 1 );
	}

	#[test]
	fn filter_build_correct() {
		let shape = Shape{ m : 60, k : 2 };
		let proto = SimpleProto::new( 1 );
		let bloomfilter  = proto.build( &shape );
		assert_eq!( bloomfilter.as_vec().len(), 2 );
		assert_eq!( bloomfilter.as_vec()[0], 3 );
		assert_eq!( bloomfilter.as_vec()[1], 0 );
		assert_eq!( bloomfilter.hamming_value(), 2);
		assert!( bloomfilter.estimate_n()-1.0 < 0.05 );
		// filter always contains itself
		assert!( bloomfilter.contains( &bloomfilter ).unwrap());
		let empty_filter : BloomFilterType = Box::new(Simple{ shape : shape.clone(), buffer : vec![0 as u32 ; 2 as usize ]}) ;
		// a filter always contains the empty filter
		assert!( bloomfilter.contains( &empty_filter ).unwrap());
		// an empty filter never contains a populated filter
		assert!( ! empty_filter.contains( &bloomfilter ).unwrap());
		// an empty filter always contains itself
		assert!( empty_filter.contains( &empty_filter ).unwrap());
	}

	#[test]
	fn shape_used_multiple_times() {
		let shape = Shape{ m : 60, k : 2 };
		let proto = SimpleProto::new( 1 );
		let bloomfilter = proto.build( &shape );
		let bloomfilter2 = proto.build( &shape );
		assert_eq!( bloomfilter.as_vec().len(), 2 );
		assert_eq!( bloomfilter.as_vec()[0], 3 );
		assert_eq!( bloomfilter.as_vec()[1], 0 );
		assert_eq!( bloomfilter2.as_vec().len(), 2 );
		assert_eq!( bloomfilter2.as_vec()[0], 3 );
		assert_eq!( bloomfilter2.as_vec()[1], 0 );
	}

	#[test]
	fn proto_collection() {
		let shape = Shape{ m : 60, k : 2 };
		// this proto will turn on the left 'k' most bits
		let proto = SimpleProto::new( 1 );
		// this proto will turn on the every other bit for a total of 'k' bits
		let proto2 = SimpleProto::new( 2 );
		let mut collection = ProtoCollection::new();
		collection.add( Box::new( proto ) );
		collection.add( Box::new( proto2 ) );
		assert_eq!( collection.len(), 2 );
		let bloomfilter = collection.build( &shape );

		assert_eq!( bloomfilter.as_vec().len(), 2 );
		assert_eq!( bloomfilter.as_vec()[0], 7 );
		assert_eq!( bloomfilter.as_vec()[1], 0 );

		//
		// test collection containing a collection
		//

		// this should yeild bits 0 and 16 (65536)
		let proto3 = SimpleProto::new( 0x100 );
		let mut collection2 = ProtoCollection::new();
		collection2.add( Box::new( collection ));
		collection2.add( Box::new( proto3 ));
		let bloomfilter2 = collection2.build( &shape );
		assert_eq!( bloomfilter2.as_vec().len(), 2 );
		assert_eq!( bloomfilter2.as_vec()[0], 7+65536 );
		assert_eq!( bloomfilter2.as_vec()[1], 0 );


	}

}
