use crate::bloomfilter::*;

/// A bloom filter that stores the bits in a BitVector
#[derive(Debug)]
pub struct Counting {
    shape: Shape,
    buffer: Vec<u64>,
}

impl Counting {
    /// Create an empty Counting instance.
    ///
    /// # param
    /// shape - The shape of the Bloom filter.
    pub fn empty_instance(shape: &Shape) -> BloomFilterType {
        let v : Vec<u64> = vec![0;shape.number_of_buckets()];
        Box::new(Counting {
            shape: shape.clone(),
            buffer: v,
        })
    }

    /// Create a Counting instance from a proto type.
    ///
    /// # param
    /// shape - The shape of the Bloom filter.
    /// proto - The prototype for the filter.
    pub fn instance(shape: &Shape, proto: &dyn Proto) -> BloomFilterType {
        let mut v : Vec<u64> =  vec![0;shape.number_of_buckets()];
        proto.bits(&shape).iter().for_each(|s| {
            v[*s] = v[*s]+1;
        });
        let simple = Counting {
            shape: shape.clone(),
            buffer: v,
        };
        return Box::new(simple);
    }

    /// Adds another Counting filter to this one.  This is not the same as a `merge_inplace`
    /// operation as
    /// the counts for each cell in the other Counting filter are added to the counts in this
    /// Counting filter.
    ///
    /// # Param
    /// other - The Counting filter to add to this one.
    pub fn add(&mut self, other : &Counting)  -> Result<bool, &str> {
        let x = &other.buffer;
        if !Counting::can_increment(&self.buffer, &x) {
            Err("Can not increment counts")
        } else {
            let enum_iter = x.iter().enumerate().filter(|(_, b)| **b > 0);
            for (a, b) in enum_iter {
                self.buffer[a] += b;
            }
            Ok(true)
        }
    }


    fn can_increment( mine : &Vec<u64>, buffer : &Vec<u64>) ->bool {
         buffer.iter().enumerate()
            .filter( |(_a,b)| **b > 0)
            .all( |(a,b) | mine[a] < (u64::MAX - b) )
    }

    fn can_decrement(mine : &Vec<u64>, buffer : &Vec<u64>) -> bool {
        buffer.iter().enumerate()
            .filter( |(_a,b)| **b > 0)
            .all( |(a,b) | mine[a] < *b )

    }

    /// Subtracts another Counting filter from this one.  This is not the same as a `remove`
    /// operation as
    /// the counts for each cell in the other Counting filter are subtracted to the counts in this
    /// Counting filter.
    ///
    /// # Param
    /// other - The Counting filter to subtract from this one.
    pub fn subtract(&mut self, other : &Counting)  -> Result<bool, &str> {
        let x = &other.buffer;
        if ! Counting::can_decrement(&self.buffer, &x) {
            Err("Can not decrement counts")
        } else {
            x.iter().enumerate().filter(|(_, b)| **b > 0).for_each(|(a, b)| self.buffer[a] -= b);
            Ok(true)
        }
    }

    fn can_increment_standard(mine : &Vec<u64>, buffer : &Vec<usize>) -> bool {
        let limit = u64::MAX - 1;
        buffer.iter().all(|idx| mine[*idx] < limit)
    }

    fn can_decrement_standard(mine : &Vec<u64>, buffer : &Vec<usize>) -> bool {
        buffer.iter().all( |idx| mine[*idx] == 0 )
    }
    /// Removes (deletes) a merged Bloom filter from this one.  This is the inverse of a
    /// `merge_inplace`  operation.
    pub fn remove(&mut self, other: &BloomFilterType) -> Result<(), &str> {
        let other_vec: &Vec<usize> = &other.indices();
        if  ! Counting::can_decrement_standard(&self.buffer, other_vec) {
            Err( "Can not decrement counters")
        } else {
            other_vec.iter().for_each(|s| self.buffer[*s] -= 1);
            Ok(())
        }
    }

}

impl Clone for Counting {
    fn clone(&self) -> Counting {
        Counting {
                 shape: self.shape().clone(),
                 buffer: self.buffer.to_vec(),
             }
    }
}

impl BloomFilter for Counting {
    /// merge this filter with another to make a new filter
    fn merge(&self, other: &BloomFilterType) -> Result<BloomFilterType, &str> {

        let mut result = self.clone();
        // Counting {
        //     shape: self.shape().clone(),
        //     buffer: self.buffer.to_vec(),
        // };
        result.merge_inplace( &other )?;
        Ok(Box::new(result) )
    }

    fn merge_inplace<'a>(&mut self, other: &BloomFilterType) -> Result<(), &'a str> {
        let other_vec: &Vec<usize> = &other.indices();
        if ! Counting::can_increment_standard(&self.buffer, other_vec) {
            Err( "Can not increment counters")
        }  else {
            other_vec.iter().for_each(|s| self.buffer[*s] += 1);
            Ok(())
        }
    }

    fn is_sparse(&self) -> bool {
        true
    }

    /// return the filter as a BitVector.
    fn vector(&self) -> Rc<BitVector> {
        let x = &self.buffer;
        let mut bv = BitVector::new( self.shape.m );
        x.iter().enumerate().filter( |(_a,b)| ** b > 0 )
            .for_each( |(a,_)| {bv.insert( a );});
        Rc::new( bv )
    }

    /// return a list of the bits that are turned on.
    fn indices(&self) -> Rc<Vec<usize>> {
        let x = &self.buffer;
        let new_vec = x.iter().enumerate().filter(|(_, b)| **b > 0)
            .map(|(a, _)| a).collect();
        Rc::new( new_vec )
    }

    fn shape(&self) -> &Shape {
        &self.shape
    }

    fn hamming_value(&self) -> usize {
        let x = &self.buffer;
        x.iter().enumerate().filter( |(_,b)| **b > 0 ).count()
    }
}
