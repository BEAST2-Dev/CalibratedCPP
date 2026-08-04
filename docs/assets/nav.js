// Single-page layout: each <section class="doc-section"> is a top-level nav item;
// its h2/h3 headings become a nested table of contents that expands when the
// section is the one currently in view. Scroll-spy highlights both.
const NAV_SECTIONS = [
  { id: "model", label: "Model" },
  { id: "constraint-tree", label: "Constraint Tree Format" },
  { id: "beast-xml", label: "Specify in BEAST3 XML" },
  { id: "beauti", label: "Specify in BEAUti" },
  { id: "lphy", label: "Specify in an LPhy script" },
  { id: "citing", label: "Citing" },
  { id: "license", label: "License"},
];

const SPY_OFFSET = 150;

function slugify(text) {
  return text
    .toLowerCase()
    .trim()
    .replace(/[^\w\s-]/g, "")
    .replace(/\s+/g, "-");
}

// Give every h2/h3 a unique id so the TOC can anchor to it (preserving any
// id already set in the HTML, e.g. specifying-calibrations).
function ensureHeadingIds() {
  const used = new Set();
  document.querySelectorAll("main h2, main h3").forEach((h) => {
    const base = h.id || slugify(h.textContent);
    let unique = base;
    let n = 1;
    while (
      used.has(unique) ||
      (!h.id && document.getElementById(unique) && document.getElementById(unique) !== h)
    ) {
      unique = `${base}-${n++}`;
    }
    h.id = unique;
    used.add(unique);
  });
}

function renderSidebar() {
  const mount = document.getElementById("sidebar-nav");
  if (!mount) return;

  ensureHeadingIds();

  const nav = document.createElement("nav");
  const list = document.createElement("ul");

  NAV_SECTIONS.forEach((section) => {
    const sectionEl = document.getElementById(section.id);
    if (!sectionEl) return;

    const li = document.createElement("li");
    li.dataset.section = section.id;

    const link = document.createElement("a");
    link.href = `#${section.id}`;
    link.className = "section-link";
    link.textContent = section.label;
    link.dataset.sectionLink = section.id;
    li.appendChild(link);

    const headings = [...sectionEl.querySelectorAll("h2, h3")];
    if (headings.length) {
      const toc = document.createElement("ul");
      toc.className = "toc";
      headings.forEach((h) => {
        const tocLi = document.createElement("li");
        tocLi.className = `toc-${h.tagName.toLowerCase()}`;
        const a = document.createElement("a");
        a.href = `#${h.id}`;
        a.textContent = h.textContent;
        a.dataset.tocId = h.id;
        tocLi.appendChild(a);
        toc.appendChild(tocLi);
      });
      li.appendChild(toc);
    }

    list.appendChild(li);
  });

  nav.appendChild(list);

  mount.innerHTML = `
    <a class="site-title" href="#model">Calibrated CPP</a>
    <span class="site-desc">A BEAST3 package for the Calibrated Coalescent Point Process</span>
  `;
  mount.appendChild(nav);
  const gh = document.createElement("a");
  gh.className = "github-link";
  gh.href = "https://github.com/BEAST2-Dev/CalibratedCPP";
  gh.target = "_blank";
  gh.rel = "noopener";
  gh.innerHTML = "View on GitHub &rarr;";
  mount.appendChild(gh);

  setupScrollSpy();
}

function setupScrollSpy() {
  const sections = NAV_SECTIONS.map((s) => document.getElementById(s.id)).filter(Boolean);
  const headings = [...document.querySelectorAll("main h2, main h3")];

  const sectionLis = new Map();
  const sectionLinks = new Map();
  document.querySelectorAll("#sidebar-nav li[data-section]").forEach((li) => {
    sectionLis.set(li.dataset.section, li);
  });
  document.querySelectorAll("#sidebar-nav a[data-section-link]").forEach((a) => {
    sectionLinks.set(a.dataset.sectionLink, a);
  });
  const headingLinks = new Map();
  document.querySelectorAll("#sidebar-nav a[data-toc-id]").forEach((a) => {
    headingLinks.set(a.dataset.tocId, a);
  });

  let ticking = false;
  function update() {
    ticking = false;
    if (!sections.length) return;

    let activeSection = sections[0].id;
    for (const s of sections) {
      if (s.getBoundingClientRect().top <= SPY_OFFSET) activeSection = s.id;
      else break;
    }

    let activeHeading = headings.length ? headings[0].id : null;
    for (const h of headings) {
      if (h.getBoundingClientRect().top <= SPY_OFFSET) activeHeading = h.id;
      else break;
    }

    // At the very bottom, force the last section/heading active.
    if (window.innerHeight + window.scrollY >= document.body.offsetHeight - 4) {
      activeSection = sections[sections.length - 1].id;
      if (headings.length) activeHeading = headings[headings.length - 1].id;
    }

    sectionLis.forEach((li, id) => li.classList.toggle("section-active", id === activeSection));
    sectionLinks.forEach((a, id) => a.classList.toggle("active", id === activeSection));
    headingLinks.forEach((a, id) => a.classList.toggle("active", id === activeHeading));
  }

  window.addEventListener(
    "scroll",
    () => {
      if (!ticking) {
        ticking = true;
        requestAnimationFrame(update);
      }
    },
    { passive: true }
  );
  window.addEventListener("resize", update);
  update();
}

function setupMobileToggle() {
  const toggle = document.getElementById("nav-toggle");
  const sidebar = document.getElementById("sidebar");
  if (!toggle || !sidebar) return;
  toggle.addEventListener("click", () => sidebar.classList.toggle("open"));
  sidebar.addEventListener("click", (e) => {
    if (e.target.closest("a")) sidebar.classList.remove("open");
  });
}

document.addEventListener("DOMContentLoaded", () => {
  renderSidebar();
  setupMobileToggle();
});
