
use std::collections::HashSet;
use crate::bloomfilter::index_producer::{IndexProducer, IndexProducerType};
use crate::bloomfilter::Shape;

pub type HasherType = Box<dyn Hasher>;

pub trait Hasher  {
    /// get the indices produced by the Hasher.
    fn indices(&self, shape: &Shape) -> Vec<IndexProducerType>;

    /// get the number of items in the Hasher
    fn size(&self) -> i32;

    /// return `true` if the the number of items in the hasher is zero.
    fn is_empty(&self) -> bool {
        self.size() == 0
    }
}


pub struct SimpleHasher {
    initial : u64,
    increment : u64,
}

impl SimpleHasher {
    pub fn new(init :u64, incr : u64) -> HasherType {
        Box::new( SimpleHasher{
            initial : init,
            increment : incr,
        } )
    }


}

struct SimpleHasherIndexProducer {
    initial : u64,
    increment : u64,
    shape : Shape,
}


impl IndexProducer for SimpleHasherIndexProducer {


    fn get_indices(&self) -> HashSet<i32> {
        let mut set:HashSet<i32> = HashSet::new();

        /*
                * Essentially this is computing a wrapped modulus from a start point and an
                * increment. So actually you only need two modulus operations before the loop.
                * This avoids any modulus operation inside the while loop. It uses a long index
                * to avoid overflow.
                */
        let mut index = (self.initial % self.shape.m as u64) as i32;
        let incr = (self.increment %  self.shape.m as u64) as i32;

        for _func_count in 0..self.shape.k {
            set.insert(index);
            index += incr;
            index = if index as usize >= self.shape.m  {
                index - self.shape.m as i32
            } else {
                index
            };
        };
        set
    }
}

impl Hasher for SimpleHasher {
    fn indices(&self, shape: &Shape) -> Vec<IndexProducerType> {
        vec![ Box::new(SimpleHasherIndexProducer {
            initial : self.initial,
            increment : self.increment,
            shape : shape.clone(),
        })]
    }

    fn size(&self) -> i32 { 1 }
}

pub struct HasherCollection {
    hashers : Vec<HasherType>,
}

impl HasherCollection {

    pub fn new() -> HasherCollection {
        HasherCollection {
            hashers: vec![],
        }
    }

    pub fn add( &mut self, hasher : HasherType ) {
        self.hashers.push( hasher );
    }
}

impl Hasher for HasherCollection {
    fn indices(&self, shape: &Shape) -> Vec<IndexProducerType> {
        let mut result : Vec<IndexProducerType> = Vec::with_capacity( self.hashers.len());
        for hasher in &self.hashers {
            for producer in hasher.indices( shape ) {
                result.push( producer );
            }
        }
        result
    }

    fn size(&self) -> i32 {
        let mut result = 0;
        for hasher in &self.hashers {
            result += hasher.size();
        };
        result
    }
}


