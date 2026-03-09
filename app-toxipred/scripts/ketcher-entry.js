// Entry point for the standalone ketcher esbuild bundle.
// Exposes React, ReactDOM, Editor (ketcher-react), and StandaloneStructServiceProvider (ketcher-standalone)
// as properties of the global __ketcher object.

import React from 'react';
import ReactDOM from 'react-dom/client';
import { Editor } from 'ketcher-react';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import 'ketcher-react/dist/index.css';

export { React, ReactDOM, Editor, StandaloneStructServiceProvider };
