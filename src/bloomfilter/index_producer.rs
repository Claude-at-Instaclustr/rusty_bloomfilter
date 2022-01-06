use std::collections::HashSet;

pub type IndexProducerType = Box<dyn IndexProducer>;

pub trait IndexProducer {
    fn get_indices(&self) -> HashSet<i32>;
}
