
pub type BitMapProducerType = Box<dyn BitMapProducer>;

pub trait BitMapProducer {
    fn get_bitmaps(&self) -> Vec<u64>;
}
