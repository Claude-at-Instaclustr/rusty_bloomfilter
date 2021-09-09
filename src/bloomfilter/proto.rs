use crate::bloomfilter::*;

/// A prototype Bloom filer.
///
/// Proto bloom filters can build filters of any shape.
/// Proto filters process the hash function to build the filter based on shape.
pub trait Proto {
    /// Creates an iterator of bits to enable based on the shape.
    ///
    /// # Arguments
    /// * `shape` - The shape of the resulting Bloom filter.
    fn bits(&self, shape: &Shape) -> Vec<usize>;

    /// Returns the number of items in this proto.
    fn size(&self) -> usize;
}

/// A Proto implementation that is a collection of Protos.  That acts as a single proto
///
/// This represents multiple objects to be placed in the filter.
pub struct ProtoCollection {
    inner: Vec<Box<dyn Proto>>,
}

/// A Proto implementation that is a single object to be placed in the filter.
#[derive(Debug)]
pub struct SimpleProto {
    start: u64,
    incr: u64,
}
/// The simple prototype
impl SimpleProto {
    /// Creates a new Simple proto from a u128
    ///
    /// # Arguments
    /// * `hash` - A 128bit hash value as a u128.
    pub fn from_u128(hash: u128) -> SimpleProto {
        let start = (hash >> 64) as u64;
        let incr = hash as u64;

        return SimpleProto { start, incr };
    }

    /// Creates a new Simple proto from a u64.
    ///
    /// # Arguments
    /// * `hash` - A 64bit hash value as a u64.
    pub fn from_u64(hash: u64) -> SimpleProto {
        let start = (hash >> 32) as u64;
        let incr = (hash & 0xFFFFFFFF) as u64;

        return SimpleProto { start, incr };
    }

    /// Creates a new Simple proto from a u32.
    ///
    /// # Arguments
    /// * `hash` - A 32bit hash value as a u32.
    pub fn from_u32(hash: u32) -> SimpleProto {
        let start = (hash >> 16) as u64;
        let incr = (hash & 0xFFFF) as u64;
        return SimpleProto { start, incr };
    }
}

impl Proto for SimpleProto {
    fn bits(&self, shape: &Shape) -> Vec<usize> {
        let mut v: Vec<usize> = Vec::with_capacity(shape.k);
        let mut accumulator: u64 = self.start;

        for _i in 0..shape.k {
            v.push((accumulator % shape.m as u64) as usize);
            accumulator += self.incr;
        }
        return v;
    }

    fn size(&self) -> usize {
        1
    }
}

impl ProtoCollection {
    /// Creates an empty collection
    pub fn new() -> ProtoCollection {
        return ProtoCollection { inner: Vec::new() };
    }

    /// Add prototypes to the collection
    ///
    /// # Arguments
    /// * `proto` - A proto to add.  May be a ProtoCollection.
    pub fn add(&mut self, proto: Box<dyn Proto>) {
        self.inner.push(proto);
    }

    /// Gets the number of `Proto`s in the collection.  This is not the same as `size()` as
    /// `size()` return the number of `SimpleProto`'s in the collection.
    pub fn count(&self) -> usize {
        return self.inner.len();
    }
}

impl Proto for ProtoCollection {
    fn bits(&self, shape: &Shape) -> Vec<usize> {
        return self.inner.iter().flat_map(|s| s.bits(shape)).collect();
    }

    fn size(&self) -> usize {
        self.inner.iter().map(|x| x.size()).sum()
    }
}

#[cfg(test)]
mod tests {

    use crate::bloomfilter::proto::*;

    #[test]
    fn simple_proto_correct() {
        let proto = SimpleProto::from_u32(1);
        assert_eq!(proto.start, 0);
        assert_eq!(proto.incr, 1);
    }

    #[test]
    fn simple_proto_bits() {
        let proto = SimpleProto::from_u32(1);
        let shape = Shape { m: 10, k: 2 };
        let v = proto.bits(&shape);

        assert_eq!(v.len(), 2);
        assert_eq!(v[0], 0);
        assert_eq!(v[1], 1);
    }

    #[test]
    fn proto_collection() {
        let shape = Shape { m: 60, k: 2 };
        // this proto will turn on the left 'k' most bits
        let proto = SimpleProto::from_u32(1);
        // this proto will turn on the every other bit for a total of 'k' bits
        let proto2 = SimpleProto::from_u32(2);
        let mut collection = ProtoCollection::new();
        collection.add(Box::new(proto));
        collection.add(Box::new(proto2));
        assert_eq!(collection.count(), 2);
        let bloomfilter = Simple::instance(&shape, &collection);
        assert_eq!(*bloomfilter.indices(), [0, 1, 2]);

        //
        // test collection containing a collection
        //

        // this should yield bits 0 and 16 (65536)
        let proto3 = SimpleProto::from_u32(0x100);
        let mut collection2 = ProtoCollection::new();
        collection2.add(Box::new(collection));
        collection2.add(Box::new(proto3));
        let bloomfilter2 = Simple::instance(&shape, &collection2);
        assert_eq!(*bloomfilter2.indices(), [0, 1, 2, 16]);
    }
}
