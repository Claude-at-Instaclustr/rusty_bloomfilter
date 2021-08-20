
/// A bloom filter that stores the bits in an array of u32 bit buffers
pub struct Simple {
	shape : Shape,
  	buffer : Vec::<u32>,	
}

impl Simple {

	/// Construct an empty Simple Bloom filter	
	pub fn new( shape : &Shape ) -> BloomFilter {
		
		let size = if shape.m % 32 > 0  {1+(shape.m / 32)} else {shape.m / 32};
		let buff = vec![0 as u32 ; size as usize ];
		BloomFilter{ shape: shape.clone(), buffer : buff }
	}
	
}

impl Bloomfilter for Simple {
 	fn asVec(&self) -> &Vec::<u32> {
		return &buffer;
	}
	
	fn shape(&self) -> &Shape {
		reeturn &shape;
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
	
	fn contains( &self, filter : &BloomFilter ) -> Result<bool, &str>
	{
		if self.shape.m != filter.shape.m || self.shape.k != filter.shape.k {
			return Err( "Shapes do not match" )
		} 
		for i in 0..self.buffer.len() {
			if self.buffer[i] & filter.buffer[i] != filter.buffer[i] {
				return Ok(false);
			}
		}
		Ok(true)
	}
	
	fn hamming_value( &self ) -> u32 {
		return self.buffer.iter().map(|i| i.count_ones() ).sum();
	}
	
	fn estimate( &self, count : u32 ) -> f32 {
		let c = count as f64;
		let m = self.shape.m as f64;
		let k = self.shape.k as f64;
		let result = -(m/k) * (1.0 - (c/m)).ln();
		return result as f32; 
	}
	
	
}