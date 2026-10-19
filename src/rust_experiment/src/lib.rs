use quick_xml::events::Event;
use quick_xml::Reader;
use std::slice;

/// Count XML elements through a small C ABI boundary.
///
/// Returns -1 for invalid input or malformed XML.
///
/// The caller must provide a readable buffer of `len` bytes for this call.
#[no_mangle]
pub unsafe extern "C" fn openms_rust_count_elements(data: *const u8, len: usize) -> isize {
    if (data.is_null() && len != 0) || len > isize::MAX as usize {
        return -1;
    }

    let bytes = if len == 0 {
        &[]
    } else {
        slice::from_raw_parts(data, len)
    };
    let mut reader = Reader::from_reader(bytes);
    let mut count = 0isize;
    let mut depth = 0usize;

    loop {
        match reader.read_event() {
            Ok(Event::Start(_)) => {
                count = match count.checked_add(1) {
                    Some(next_count) => next_count,
                    None => return -1,
                };
                depth = match depth.checked_add(1) {
                    Some(next_depth) => next_depth,
                    None => return -1,
                };
            }
            Ok(Event::Empty(_)) => {
                count = match count.checked_add(1) {
                    Some(next_count) => next_count,
                    None => return -1,
                };
            }
            Ok(Event::End(_)) => {
                if depth == 0 {
                    return -1;
                }
                depth -= 1;
            }
            Ok(Event::Eof) => return if depth == 0 { count } else { -1 },
            Ok(_) => {}
            Err(_) => return -1,
        }
    }
}
