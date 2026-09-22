(function () {
  // Local copies first so the deck presents without a network; CDNs remain as a fallback.
  const CDNS = ['vendor/rive-canvas.js', 'https://unpkg.com/@rive-app/canvas@2.21.6', 'https://cdn.jsdelivr.net/npm/@rive-app/canvas@2.21.6'];
  const LOCAL_WASM = 'vendor/rive.wasm';
  let loading = null;
  function loadRuntime() {
    if (window.rive) return Promise.resolve();
    if (loading) return loading;
    loading = new Promise((resolve, reject) => {
      const tryIdx = (i) => {
        if (i >= CDNS.length) return reject(new Error('runtime unavailable'));
        const s = document.createElement('script'); s.src = CDNS[i]; s.async = true;
        const t = setTimeout(() => { s.remove(); tryIdx(i + 1); }, 8000);
        s.onload = () => {
          clearTimeout(t);
          try { window.rive.RuntimeLoader.setWasmUrl(new URL(LOCAL_WASM, location.href).href); } catch (e) {}
          resolve();
        };
        s.onerror = () => { clearTimeout(t); s.remove(); tryIdx(i + 1); };
        document.head.appendChild(s);
      };
      tryIdx(0);
    });
    return loading;
  }
  class RivePlayer extends HTMLElement {
    connectedCallback() {
      if (this._init) return; this._init = true;
      this.style.display = 'block'; this.style.width = '100%'; this.style.height = '100%'; this.style.cursor = 'pointer'; this.style.position = 'relative';
      const c = this._canvas = document.createElement('canvas');
      c.style.width = '100%'; c.style.height = '100%'; c.style.display = 'block';
      this.appendChild(c);
      const msg = document.createElement('div');
      msg.style.cssText = 'position:absolute;left:0;top:0;font:14px ui-monospace,monospace;color:#736F7D;padding:8px';
      msg.textContent = 'loading animation runtime…';
      this.appendChild(msg);
      // steps="Setup,Finish B,Finish C,Finish A" → click advances; each step plays once and holds its last frame
      this._steps = (this.getAttribute('steps') || '').split(',').map(s => s.trim()).filter(Boolean);
      // loops="Running ABC,Running AC,Running A," → looping animation layered on after step N finishes (blank = none)
      this._loops = (this.getAttribute('loops') || '').split(',').map(s => s.trim());
      this._i = -1;
      const start = () => {
        const opts = { src: this.getAttribute('src'), canvas: c, autoplay: false,
          artboard: this.getAttribute('artboard') || undefined,
          layout: new rive.Layout({ fit: rive.Fit.Contain, alignment: rive.Alignment.Center }),
          onLoad: () => {
            this._fit(); msg.remove(); this._loaded = true;
            try { this._log('animations: ' + r.animationNames.join(' | ') + '  state machines: ' + r.stateMachineNames.join(' | ')); } catch (e) {}
            if (this._isActiveSlide()) this._begin();
          },
          onStop: (ev) => {
            // a one-shot step finished → start its idle loop, if any
            const names = (ev && ev.data) || [];
            const cur = this._steps[this._i], loop = this._loops[this._i];
            if (this._i >= 0 && loop && names.indexOf(cur) !== -1 && this._rive) { this._rive.play(loop); this._log('loop: ' + loop); }
          },
          onLoadError: (e) => { msg.textContent = 'could not load .riv'; this._log('load error ' + e); } };
        const sm = this.getAttribute('state-machine') || this.getAttribute('statemachine');
        const an = this.getAttribute('animation');
        if (this._steps.length) opts.animations = this._steps[0];
        else if (sm) opts.stateMachines = sm;
        else if (an) opts.animations = an;
        const r = this._rive = new rive.Rive(opts);
        let rt; const refit = () => { clearTimeout(rt); rt = setTimeout(() => this._fit(), 100); };
        new ResizeObserver(refit).observe(c);
        window.addEventListener('resize', refit);
        const stage0 = this.closest('deck-stage'); if (stage0) stage0.addEventListener('slidechange', refit);
        this.addEventListener('click', () => this._steps.length ? this.next() : this.reset());
        this._key = (e) => {
          if (!this._loaded || !this._steps.length || !this._isActiveSlide()) return;
          if (e.metaKey || e.ctrlKey || e.altKey) return;
          const t = e.composedPath ? e.composedPath()[0] : e.target;
          if (t && (t.isContentEditable || /^(INPUT|TEXTAREA|SELECT)$/.test(t.tagName))) return;
          // swallow → / space while steps remain, ← while steps have played; otherwise let the deck navigate
          if ((e.key === 'ArrowRight' || e.key === ' ') && this._i < this._steps.length - 1) { this.next(); }
          else if (e.key === 'ArrowLeft' && this._i >= 0) { this.back(); }
          else if (e.key === 'r') { this.restart(); }
          else return;
          e.preventDefault(); e.stopImmediatePropagation();
        };
        window.addEventListener('keydown', this._key, true);
        // inside a deck: (re)start the sequence whenever this slide becomes active
        const stage = this.closest('deck-stage');
        if (stage) stage.addEventListener('slidechange', (e) => {
          const active = e.detail && e.detail.slide ? e.detail.slide.contains(this) : this._isActiveSlide();
          if (active && this._loaded) { this._started ? this.restart() : this._begin(); }
        });
      };
      const go = () => loadRuntime().then(start, () => { msg.textContent = 'animation runtime unavailable (no network to CDN)'; });
      if (document.readyState === 'complete') setTimeout(go, 0); else window.addEventListener('load', go, { once: true });
    }
    // deck-stage scales slides with a CSS transform, so the canvas's CSS size understates its on-screen size;
    // fold the transform scale into the device pixel ratio so the drawing surface matches real pixels
    _fit() {
      const r = this._rive, c = this._canvas; if (!r || !c) return;
      const rect = c.getBoundingClientRect();
      const scale = c.clientWidth ? rect.width / c.clientWidth : 1;
      r.resizeDrawingSurfaceToCanvas((window.devicePixelRatio || 1) * scale);
    }
    _isActiveSlide() {
      const stage = this.closest('deck-stage');
      if (!stage) return true;
      const sec = this.closest('section');
      if (!sec) return true;
      const cs = getComputedStyle(sec);
      return cs.display !== 'none' && cs.visibility !== 'hidden' && sec.getAttribute('aria-hidden') !== 'true';
    }
    _begin() {
      const r = this._rive; this._started = true;
      if (!this._steps.length) { r.play(); return; }
      // show the first frame of the first step and wait for a click
      this._i = -1;
      try { r.stop(); } catch (e) {}
      r.play(this._steps[0]); r.pause(this._steps[0]);
      this._log('ready — click to start');
    }
    _playStep() {
      const r = this._rive, name = this._steps[this._i];
      // stop whatever is running; play this step once
      try { r.stop(); } catch (e) {}
      r.play(name);
      this._log('step ' + (this._i + 1) + '/' + this._steps.length + ': ' + name);
    }
    next() {
      if (!this._rive || !this._steps.length) return;
      if (this._i < this._steps.length - 1) { this._i++; this._playStep(); }
    }
    back() {
      if (!this._rive || this._i < 0) return;
      // undo the last step: show the first frame of the step just undone
      const r = this._rive, name = this._steps[this._i]; this._i--;
      try { r.stop(); } catch (e) {}
      r.play(name); r.pause(name);
      this._log('back to before: ' + name);
    }
    restart() { if (!this._rive) return; this._rive.reset({ autoplay: false, artboard: this.getAttribute('artboard') || undefined, animations: this._steps[0] }); setTimeout(() => this._begin(), 0); }
    _log(m) { console.log('[rive-player] ' + m); }
    play() { this._rive && this._rive.play(); }
    reset() { if (!this._rive) return; const an = this.getAttribute('animation'); this._rive.reset({ autoplay: true, animations: an || undefined }); }
  }
  customElements.define('rive-player', RivePlayer);
})();
