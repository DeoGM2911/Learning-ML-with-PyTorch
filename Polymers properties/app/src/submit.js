// Navbar toggle (same as index.html)
    const toggle = document.getElementById('menu-toggle');
    const navLinks = document.getElementById('nav-links');
    toggle.addEventListener('click', () => {
      navLinks.classList.toggle('open');
      toggle.classList.toggle('active');
    });

const chatForm  = document.getElementById('chat-form');
const textInput = document.getElementById('text-input');
const fileInput = document.getElementById('file-input');
const chatWindow= document.getElementById('chat-window');

// 1) ‘Enter’ to send
textInput.addEventListener('keydown', e => {
  if (e.key === 'Enter' && !e.shiftKey) {
    e.preventDefault();
    chatForm.requestSubmit();
  }
});

// 2) Auto-submit on file select
fileInput.addEventListener('change', () => {
  if (fileInput.files.length > 0) {
    // Optional: clear the text field so only the image sends
    textInput.value = '';
    chatForm.requestSubmit();
  }
});


chatForm.addEventListener('submit', async e => {
  e.preventDefault();
  const text = textInput.value.trim();
  const file = fileInput.files[0];
  if (!text && !file) return;

  // 1) render the “pending” user bubble
  const userMsg = document.createElement('div');
  userMsg.className = 'message user';
  userMsg.textContent = text || file.name;
  chatWindow.appendChild(userMsg);
  chatWindow.scrollTop = chatWindow.scrollHeight;

  // 2) send to /api/predict
  const formData = new FormData();
  if (file) formData.append('file', file);
  if (text) formData.append('text', text);
  textInput.value = '';
  fileInput.value = '';

  let taskId;
  try {
    const resp = await fetch('/api/predict', {
      method: 'POST',
      body: formData
    });
    const { task_id, status } = await resp.json();
    if (status !== 'pending') throw new Error('Unexpected status');
    taskId = task_id;
  } catch (err) {
    console.error('Error enqueuing task:', err);
    return;
  }

  // 3) poll /api/result/<taskId>
  async function pollResult() {
    try {
      const res = await fetch(`/api/result/${taskId}`);
      const payload = await res.json();

      if (payload.status === 'pending') {
        // not ready yet → try again in 1s
        setTimeout(pollResult, 1000);
      } else if (payload.status === 'done') {
        // 4) render the bot reply
        const botMsg = document.createElement('div');
        botMsg.className = 'message bot';
        botMsg.textContent = payload.reply;
        chatWindow.appendChild(botMsg);
        chatWindow.scrollTop = chatWindow.scrollHeight;
      } else {
        console.error('Task failed:', payload);
      }
    } catch (err) {
      console.error('Error polling result:', err);
    }
  }

  pollResult();
});
