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
