use std::io;

/// Whether an error chain contains a broken-pipe error.  A downstream `head`
/// closing the pipe should terminate us quietly rather than print an error.
pub fn is_broken_pipe(err: &anyhow::Error) -> bool {
    err.chain().any(|cause| {
        cause
            .downcast_ref::<io::Error>()
            .is_some_and(|e| e.kind() == io::ErrorKind::BrokenPipe)
    })
}
