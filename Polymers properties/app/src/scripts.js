// scripts.js
document.addEventListener('DOMContentLoaded', () => {
  const toggle   = document.getElementById('menu-toggle');
  const navLinks = document.getElementById('nav-links');
  const nav      = document.querySelector('nav');

  toggle.addEventListener('click', e => {
    e.stopPropagation();
    navLinks.classList.toggle('open');
    toggle.classList.toggle('active');
  });

  document.addEventListener('click', e => {
    if (navLinks.classList.contains('open') && !nav.contains(e.target)) {
      navLinks.classList.remove('open');
      toggle.classList.remove('active');
    }
  });
});

// optional: standalone function can go here too
function scrollToPredictions() {
  window.location.href = '/predictions';
}

*,
*::before,
*::after {
  box-sizing: border-box;
}
html, body {
  margin: 0;
  padding: 0;
  height: 100%;
}

/* 2) Make the body a column flex container */
body {
  display: flex;
  flex-direction: column;
  /* optional: prevent horizontal overflow if something is slightly too wide */
  overflow-x: hidden;
}

/* 3) Nav and footer should take only their natural height */
nav,
footer {
  flex: none;
}

/* 4) Let your chat container fill the remaining space */
.chat-container {
  display: flex;
  flex-direction: column;
  flex: 1;               /* grow to fill */
  overflow: hidden;      /* contain the inner scroll */
}

/* 5) The inner window scrolls, not the whole page */
.chat-window {
  flex: 1;
  overflow-y: auto;
}
