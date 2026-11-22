// Runs the script of a CanSNPer2 result page under node with a minimal stand-in for the
// browser document, so the tests can see what the page draws without a browser.
// usage: node page_harness.js page.html [javascript run after the page script]
// prints JSON: {errors, html: {id: innerHTML}, rows: [row html], viewBox, downloads: [{name, type, text}]}
"use strict";
const fs = require("fs");
const vm = require("vm");
const page = fs.readFileSync(process.argv[2], "utf8");
const after = process.argv[3] || "";
const scripts = [...page.matchAll(/<script[^>]*>([\s\S]*?)<\/script>/g)].map(m => m[1]);
const ids = new Set([...page.split("<script")[0].matchAll(/id="([^"]+)"/g)].map(m => m[1]));
const errors = [];
const downloads = [];
const blobs = {};
let blobCount = 0;

class El {
	constructor(id, tag) {
		this.id = id; this.tagName = (tag || "div").toUpperCase(); this._html = ""; this.children = [];
		this.style = {}; this.attrs = {}; this.className = ""; this.value = ""; this.checked = false;
		this.listeners = {}; this.parentNode = null; this._text = "";
	}
	get textContent() { return this._text; }
	set textContent(v) { this._text = String(v); this._html = this._text.replace(/&/g, "&amp;").replace(/</g, "&lt;"); }
	get innerHTML() { return this._html; }
	set innerHTML(v) { this._html = String(v); this.children = []; this._svg = null; }
	insertAdjacentHTML(pos, v) { this._html += v; }
	appendChild(c) { this.children.push(c); c.parentNode = this; return c; }
	removeChild(c) { this.children = this.children.filter(x => x !== c); }
	insertBefore(c) { errors.push(c.textContent); this.children.unshift(c); return c; }
	setAttribute(k, v) { this.attrs[k] = String(v); }
	getAttribute(k) { return k in this.attrs ? this.attrs[k] : null; }
	hasAttribute(k) { return k in this.attrs; }
	addEventListener(t, f) { (this.listeners[t] = this.listeners[t] || []).push(f); }
	getBoundingClientRect() { return {left: 0, top: 0, width: 1000, height: 640}; }
	click() {
		if (this.onclick) { this.onclick({target: this}); }
		if (this.tagName === "A" && this.download) {
			const blob = blobs[this.href] || {text: "", type: ""};
			downloads.push({name: this.download, type: blob.type, text: blob.text});
		}
	}
}
const elements = {};
function byId(id) {
	if (!ids.has(id)) { return null; }
	if (!elements[id]) { elements[id] = new El(id); }
	return elements[id];
}
const body = new El("body", "body");
const document = {
	body: body,
	title: "",
	getElementById: byId,
	createElement: tag => new El(null, tag),
	querySelector(sel) {
		let m = sel.match(/^#(\S+) svg$/);
		if (m) {
			const el = byId(m[1]);
			if (!el || el.innerHTML.indexOf("<svg") === -1) { return null; }
			if (!el._svg) { el._svg = new El(null, "svg"); }
			return el._svg;
		}
		m = sel.match(/^#(\S+) tr\[data-i='(\d+)'\]$/);
		if (m) { return (byId(m[1]) || {children: []}).children.find(c => c.attrs["data-i"] === m[2]) || null; }
		throw new Error("page_harness.js does not know the selector " + sel);
	},
	querySelectorAll(sel) {
		const m = sel.match(/^#(\S+) tr$/);
		if (m) { return (byId(m[1]) || {children: []}).children; }
		throw new Error("page_harness.js does not know the selector " + sel);
	},
};
class Blob {
	constructor(parts, opts) { this.text = parts.join(""); this.type = (opts || {}).type || ""; }
}
const URL = {
	createObjectURL(b) { const key = "blob:" + (++blobCount); blobs[key] = b; return key; },
	revokeObjectURL() {},
};
const timers = [];
const window = {listeners: {}, addEventListener(t, f) { (this.listeners[t] = this.listeners[t] || []).push(f); }};
const context = vm.createContext({document, window, Blob, URL, console, setTimeout: f => timers.push(f), Math, JSON, String, Object, RegExp, encodeURIComponent});
context.window = window;
try {
	vm.runInContext(scripts.join("\n") + "\n" + after, context, {filename: "page.js"});
} catch (e) {
	if (window.onerror) { window.onerror(String(e.message), "", 0); } else { errors.push(String(e.stack)); }
	errors.push(String(e.stack));
}
const out = {errors, html: {}, rows: [], viewBox: null, downloads};
for (const id in elements) { out.html[id] = elements[id].innerHTML; }
if (elements.rows) { out.rows = elements.rows.children.map(c => ({cls: c.className, html: c.innerHTML})); }
if (elements.tree && elements.tree._svg) { out.viewBox = elements.tree._svg.attrs.viewBox || null; }
process.stdout.write(JSON.stringify(out));
