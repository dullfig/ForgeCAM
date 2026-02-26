use std::marker::PhantomData;

/// Type-safe index into an Arena<T>.
#[derive(Debug)]
pub struct Idx<T> {
    raw: u32,
    _marker: PhantomData<T>,
}

// Manual impls to avoid requiring T: Clone/Copy/etc.
impl<T> Clone for Idx<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Idx<T> {}
impl<T> PartialEq for Idx<T> {
    fn eq(&self, other: &Self) -> bool {
        self.raw == other.raw
    }
}
impl<T> Eq for Idx<T> {}
impl<T> std::hash::Hash for Idx<T> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.raw.hash(state);
    }
}

impl<T> Idx<T> {
    pub fn from_raw(raw: u32) -> Self {
        Self {
            raw,
            _marker: PhantomData,
        }
    }
    pub fn raw(self) -> u32 {
        self.raw
    }
}

/// Contiguous arena storing entities of type T, indexed by Idx<T>.
pub struct Arena<T> {
    items: Vec<T>,
}

impl<T> Arena<T> {
    pub fn new() -> Self {
        Self { items: Vec::new() }
    }

    pub fn alloc(&mut self, item: T) -> Idx<T> {
        let idx = Idx::from_raw(self.items.len() as u32);
        self.items.push(item);
        idx
    }

    pub fn get(&self, idx: Idx<T>) -> &T {
        &self.items[idx.raw as usize]
    }

    pub fn get_mut(&mut self, idx: Idx<T>) -> &mut T {
        &mut self.items[idx.raw as usize]
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (Idx<T>, &T)> {
        self.items
            .iter()
            .enumerate()
            .map(|(i, item)| (Idx::from_raw(i as u32), item))
    }
}

impl<T> Default for Arena<T> {
    fn default() -> Self {
        Self::new()
    }
}
