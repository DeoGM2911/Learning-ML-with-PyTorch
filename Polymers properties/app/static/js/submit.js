/* submit.js — ES module */
const { chatForm, textInput, fileInput, chatWindow } = window.__chatHandles;

/* ---------- helpers ---------- */
const escapeHTML = (str = '') =>
  str.replace(/[&<>"']/g, ch =>
    ({ '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;' }[ch]));

function addMessage({ html, sender }) {
  const div = document.createElement('div');
  div.className = `message ${sender}`;   // 'user' | 'bot'
  div.innerHTML = html;
  chatWindow.appendChild(div);
  chatWindow.scrollTop = chatWindow.scrollHeight;
}

/* ---------- UX niceties ---------- */
textInput.addEventListener('keydown', e => {
  if (e.key === 'Enter' && !e.shiftKey) {
    e.preventDefault();
    chatForm.requestSubmit();
  }
});

fileInput.addEventListener('change', () => {
  if (fileInput.files.length) {
    textInput.value = '';                // send file only
    chatForm.requestSubmit();
  }
});

/* ---------- main submit handler ---------- */
chatForm.addEventListener('submit', async e => {
  e.preventDefault();

  const text = textInput.value.trim();
  const file = fileInput.files[0];
  if (!text && !file) return;

  /* 1. optimistic user bubble ---------------------------------------- */
  let userHtml = '';
  if (file) {
    if (file.type.startsWith('image/')) {
      // show the image itself
      const url = URL.createObjectURL(file);
      userHtml = `<img src="${url}" alt="${escapeHTML(file.name)}"
                class="uploaded-img" style="max-width: 100%; height: auto; object-fit: contain;"
                onload="URL.revokeObjectURL(this.src)">`;
    } else {
      // non-image file → just the name
      userHtml = `<i class="filename">${escapeHTML(file.name)}</i>`;
    }
  } else {
    userHtml = escapeHTML(text);
  }
  addMessage({ html: userHtml, sender: 'user' });

  /* 2. send to backend ------------------------------------------------ */
  const formData = new FormData();
  if (file) formData.append('file', file);
  if (text) formData.append('text', text);
  textInput.value = '';
  fileInput.value = '';

  let taskId;
  try {
    const { task_id, status } = await fetch('/api/predict', {
      method: 'POST',
      body: formData
    }).then(r => r.json());

    if (status !== 'pending') throw new Error('enqueue failed');
    taskId = task_id;
  } catch (err) {
    console.error(err);
    addMessage({
      html: '<span style="color:red">❌ could not send message. Please check your connection!</span>',
      sender: 'bot'
    });
    return;
  }

  /* 3. show loading indicator ---------------------------------------- */
const loadingDiv = document.createElement('div');
loadingDiv.className = 'message bot loading';
loadingDiv.innerHTML = `<span class="loading-dots"><span>.</span><span>.</span><span>.</span></span>`;
chatWindow.appendChild(loadingDiv);
chatWindow.scrollTop = chatWindow.scrollHeight;

/* 4. poll for result & replace loading bubble ---------------------- */
(async function poll() {
  try {
    const { status, reply_html } =
      await fetch(`/api/result/${taskId}`).then(r => r.json());

    if (status === 'pending') {
      setTimeout(poll, 1000);
    } else if (status === 'done') {
      loadingDiv.outerHTML = '';  // remove loading message
      addMessage({ html: reply_html, sender: 'bot' });
    } else {
      throw new Error('task failed');
    }
  } catch (err) {
    console.error(err);
    loadingDiv.outerHTML = '';  // remove loading message
    addMessage({
      html: '<span style="color:red">❌ error fetching reply</span>',
      sender: 'bot'
    });
  }
})();
})