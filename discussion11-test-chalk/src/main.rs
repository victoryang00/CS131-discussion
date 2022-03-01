using std::array::IntoIter;

trait Iterator {
    type Item;
}

impl<A> Iterator for IntoIter<A> {
    type Item = A;
}

impl<T: Iterator> Iterator for Enumerate<T> {
    type Item = (usize, T::Item);
}

fn main() {
    println!("Hello, world!");
}
