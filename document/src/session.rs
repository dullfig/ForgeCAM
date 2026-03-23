use crate::Document;

/// Identifies a document within a session. Stable for the lifetime of the session.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DocumentId(u32);

/// The state of a document's background computation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ComputeStatus {
    /// No background work running.
    Idle,
    /// Feature recognition in progress.
    Recognizing,
    /// Toolpath generation in progress.
    GeneratingToolpaths,
    /// Background work finished, results ready.
    Ready,
    /// Background work failed.
    Failed(String),
}

/// A document slot in the session — the document plus its status.
#[derive(Debug)]
struct DocumentSlot {
    document: Document,
    status: ComputeStatus,
}

/// Holds all open documents. One is active (edited by the user, rendered by GUI).
/// Others may be idle or running background computation.
///
/// No cross-document references. Documents are independent.
#[derive(Debug)]
pub struct Session {
    slots: Vec<DocumentSlot>,
    active: Option<usize>,
    next_id: u32,
    ids: Vec<DocumentId>,
}

impl Session {
    pub fn new() -> Self {
        Self {
            slots: Vec::new(),
            active: None,
            next_id: 0,
            ids: Vec::new(),
        }
    }

    /// Add a document to the session and make it active. Returns its ID.
    pub fn add(&mut self, document: Document) -> DocumentId {
        let id = DocumentId(self.next_id);
        self.next_id += 1;
        let idx = self.slots.len();
        self.slots.push(DocumentSlot {
            document,
            status: ComputeStatus::Idle,
        });
        self.ids.push(id);
        self.active = Some(idx);
        id
    }

    /// Switch the active document. Returns false if the ID is not found.
    pub fn set_active(&mut self, id: DocumentId) -> bool {
        if let Some(idx) = self.index_of(id) {
            self.active = Some(idx);
            true
        } else {
            false
        }
    }

    /// Close a document and remove it from the session.
    /// If it was active, active becomes None.
    pub fn close(&mut self, id: DocumentId) -> Option<Document> {
        let idx = self.index_of(id)?;
        let slot = self.slots.remove(idx);
        self.ids.remove(idx);
        match self.active {
            Some(a) if a == idx => self.active = None,
            Some(a) if a > idx => self.active = Some(a - 1),
            _ => {}
        }
        Some(slot.document)
    }

    /// Get a reference to the active document.
    pub fn active(&self) -> Option<&Document> {
        self.active.map(|idx| &self.slots[idx].document)
    }

    /// Get a mutable reference to the active document.
    pub fn active_mut(&mut self) -> Option<&mut Document> {
        self.active.map(|idx| &mut self.slots[idx].document)
    }

    /// Get a reference to a document by ID.
    pub fn get(&self, id: DocumentId) -> Option<&Document> {
        self.index_of(id).map(|idx| &self.slots[idx].document)
    }

    /// Get the compute status of a document.
    pub fn status(&self, id: DocumentId) -> Option<&ComputeStatus> {
        self.index_of(id).map(|idx| &self.slots[idx].status)
    }

    /// Set the compute status of a document.
    pub fn set_status(&mut self, id: DocumentId, status: ComputeStatus) -> bool {
        if let Some(idx) = self.index_of(id) {
            self.slots[idx].status = status;
            true
        } else {
            false
        }
    }

    /// List all open documents: (id, name, status).
    pub fn list(&self) -> Vec<(DocumentId, &str, &ComputeStatus)> {
        self.ids
            .iter()
            .zip(self.slots.iter())
            .map(|(id, slot)| (*id, slot.document.metadata.name.as_str(), &slot.status))
            .collect()
    }

    /// Number of open documents.
    pub fn count(&self) -> usize {
        self.slots.len()
    }

    /// The active document's ID, if any.
    pub fn active_id(&self) -> Option<DocumentId> {
        self.active.map(|idx| self.ids[idx])
    }

    fn index_of(&self, id: DocumentId) -> Option<usize> {
        self.ids.iter().position(|i| *i == id)
    }
}

impl Default for Session {
    fn default() -> Self {
        Self::new()
    }
}
