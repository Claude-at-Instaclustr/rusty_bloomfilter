pub struct BitMap {}

impl BitMap {
    const DIVIDE_BY_64: i32 = 6;

    /// if the `bitmaps` array has the `idx` bit enabled returns `true`
    pub fn contains(bitmaps: &Vec<u64>, idx: &i32) -> bool {
        *idx >= 0
            && BitMap::get_long_index(idx) < bitmaps.len()
            && (bitmaps[BitMap::get_long_index(idx)] & BitMap::get_long_bit(idx)) != 0
    }

    /// returns the bit map pattern for the position within a u64.
    /// the proper bit will be enabled all other bits will be zero.
    pub fn get_long_bit(bit_index: &i32) -> u64 {
        let result :u64 = 1;
        result << (bit_index % 64)
    }

    /// returns the index of the bitmap that contains the specified bit
    pub fn get_long_index(bit_index: &i32) -> usize {
        (bit_index >> BitMap::DIVIDE_BY_64) as usize
    }

    /// returns the number of bitmaps necessary to contain `number_of_bits` bits.
    pub fn number_of_bitmaps(number_of_bits: &i32) -> usize {
        match *number_of_bits == 0 {
            true => 0,
            false => (((number_of_bits - 1) >> BitMap::DIVIDE_BY_64) + 1) as usize
        }
    }

    /// sets (enables) the bit specified by `idx` in the bitmap.
    pub fn set(bitmaps: &mut Vec<u64>, idx: &i32) {
        let block = BitMap::get_long_index(idx);
        while bitmaps.len() <= block {
            bitmaps.push( 0 );
        }
        bitmaps[block] |= BitMap::get_long_bit(idx);
    }
}
