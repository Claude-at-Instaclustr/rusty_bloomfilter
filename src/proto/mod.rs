use std::convert::TryInto;
use std::vec::Vec;
    
pub trait Proto {
	fn build( &self, shape : &Shape ) -> BloomFilter;
	
	fn add_to( &self,  filter : &mut BloomFilter );
}

struct ProtoCollection {
	inner: Vec<ProtoImpl>,
}

pub struct ProtoImpl {
	start : u32,
	incr : u32,
}

#[derive(Copy, Clone)]
pub struct Shape {
	// expected number of items
	n : u32,
	// size of array in bits
	m : u32,
	// number of functions
	k : u32,
}

pub struct BloomFilter {
	shape : Shape,
  	buffer : Vec::<u32>,	
}

impl ProtoImpl {
	pub fn new( full : u64 ) -> ProtoImpl {
		let start =   (full >> 32) as u32 ;
		let incr = full  as u32;
		
		return ProtoImpl{ start, incr };
	}
//}
//
//impl Proto for ProtoImpl {
	fn build( &self, shape : &Shape ) -> BloomFilter {
		let mut filter : BloomFilter =  BloomFilter::new( shape );
		self.add_to( &mut filter );
		return filter;
	}
	
	fn add_to( &self,  filter : &mut BloomFilter ) {
		let mut accumulator = self.start;
		for _i in 0..filter.shape.k {
			filter.on( accumulator % filter.shape.m );
			accumulator = accumulator + self.incr;
		}
	}
}

impl ProtoCollection {
	pub fn new() -> ProtoCollection {
		return ProtoCollection{ inner : Vec::new() }
	}
	
	pub fn add(&mut self,  proto : ProtoImpl ) {
		self.inner.push( proto );
	}
	
	pub fn add_collection(&mut self,  proto : ProtoCollection ) {
		for i in proto.inner {
		    self.inner.push( i );
		}
	}
	
	pub fn len(&self) -> usize {
		return self.inner.len();
	}
}

impl Proto for ProtoCollection {
	fn build( &self, shape : &Shape ) -> BloomFilter {
		let mut filter : BloomFilter =  BloomFilter::new( shape );
		self.add_to( &mut filter );
		return filter
	}
	
	fn add_to( &self, filter : &mut BloomFilter) {
		for v in &self.inner {
			v.add_to(  filter )
		}
	}
}

impl Shape {
	pub fn false_positives( &self ) -> f64 {
		let k = self.k as f64;
		let n = self.n as f64;
		let m = self.m as f64;
		return ( 1.0-((-k*n)/m).exp()).powf(k);
	}
}

impl BloomFilter {
	pub fn new( shape : &Shape ) -> BloomFilter {
		
		let size = if shape.m % 32 > 0  {1+(shape.m / 32)} else {shape.m / 32};
		let buff = vec![0 as u32 ; size as usize ];
		BloomFilter{ shape: shape.clone(), buffer : buff }
	}
	
	
	pub fn on( &mut self, bit : u32 ) {
		if bit >= self.shape.m {
			panic!( "Bit position too large")
		}
		let offset = bit % 32;
		let position : usize = (bit / 32).try_into().unwrap();
		self.buffer[position] |= 1<<offset;
	}
	
	pub fn off( &mut self, bit : u32 ) {
		if bit >= self.shape.m {
			panic!( "Bit position too large")
		}
		let offset = bit % 32;
		let position : usize = (bit / 32).try_into().unwrap() ;
		self.buffer[position] &= !(1<<offset);
	}
	
	pub fn contains( &self, filter : &BloomFilter ) -> Result<bool, &str>
	{
		if self.shape.m == filter.shape.m && self.shape.k == filter.shape.k {
			return Err( "Shapes do not match" )
		} 
		for i in 0..self.buffer.len() {
			if self.buffer[i] & filter.buffer[i] != filter.buffer[i] {
				return Ok(false);
			}
		}
		Ok(true)
	}
	
	pub fn hamming_value( &self ) -> u32 {
		return self.buffer.iter().map(|i| i.count_ones() ).sum();
	}
	
	fn estimate( &self, count : u32 ) -> f32 {
		let c = count as f64;
		let m = self.shape.m as f64;
		let k = self.shape.k as f64;
		let result = -(m/k) * (1.0 - (c/m)).ln();
		return result as f32; 
	}
	
	pub fn estimate_n( &self ) -> f32 {
		return self.estimate( self.hamming_value() );
	}
	
	pub fn estimate_union( &self,  other : &BloomFilter ) -> f32 {
		let mut count = 0;
		for i in 0..self.buffer.len() {
			let x = self.buffer[i] | other.buffer[i];
			count += x.count_ones();
		}
		return self.estimate( count );
	}

	
	pub fn estimate_intersection( &self, other : &BloomFilter ) -> f32 {
		return self.estimate_n() + other.estimate_n() - self.estimate_union( other );
	}
}

#[cfg(test)]
mod tests {
    #[test]
    fn shape_false_positives() {
		let shape : super::Shape = super::Shape{ m : 134_191, n : 4000 , k : 23};
		
        assert!(shape.false_positives()-(1.0/9_994_297.0) < 0.0000001 );
    }

	#[test]
	fn filter_protoImpl_correct() {
		let proto : super::ProtoImpl = super::ProtoImpl::new( 1 );
		assert_eq!( proto.start, 0 );
		assert_eq!( proto.incr, 1 );
	}

	#[test]
	fn filter_build_correct() {
		let shape : super::Shape = super::Shape{ m : 60, n : 4, k : 2 };
		let proto : super::ProtoImpl = super::ProtoImpl::new( 1 );
		let bloomfilter : super::BloomFilter = proto.build( &shape );
		assert_eq!( bloomfilter.buffer.len(), 2 );
		assert_eq!( bloomfilter.buffer[0], 3 );
		assert_eq!( bloomfilter.buffer[1], 0 );
	}

}