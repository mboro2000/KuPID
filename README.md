<a id="readme-top"></a>

<br />
<div align="left">
  <h1 align="left">KuPID: Kmer-based Upstream Preprocessing for Isoform Discovery</h1>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#introduction">Introduction</a> </li>
    <li><a href="#install">Install</a></li>
    <li><a href="#quickstart">Quickstart</a></li>
    <li><a href="#output">Output</a></li>
    <li><a href="#citation">Citation</a></li>
  </ol>
</details>

<!-- intro -->
## Introduction

Eukaryotic genes can encode multiple protein isoforms based on alternative splicing of exonic and intronic regions. More than 95\% of human genes have been found to undergo alternative splicing. Alternative isoforms have been implicated in a wide-range of biological processes, such cell differentiation, stress response, and tumorigenesis. Most modern novel isoform discovery methods function by identifying and assembling exon splice junctions from an RNAseq sample. However, splice junctions can only be accurately annotated with time-intensive dynamic programming alignment.

KuPID is a preprocessing method designed for RNAseq analysis of long transcript reads. When given an RNAseq sample, KuPID will filter out the reads expressed from reference (annotated) isoforms. The remaining reads can be submitted to the downstream alignment and isoform discovery software of the user's choice. KuPID operates by applying kmer sketching methods to quickly pseudo-align RNAseq reads to a reference transcriptome. The KuPID-processed samples require less time for downstream alignment, and routinely discover more true novel isoforms.

<!-- install -->
## Install

KuPID is currently available as a rust crate. 

Requirements:
1. rust programming language and associated tools such as cargo are required and assumed to be in PATH.
2. A c compiler (e.g. GCC)
3. make

```
git clone https://github.com/bluenote-1577/skani
cd skani

# If default rust install directory is ~/.cargo
cargo install --path . --root ~/.cargo

# If ~/.cargo doesn't exist use below commands instead
#cargo build --release
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>


### Built With

This section should list any major frameworks/libraries used to bootstrap your project. Leave any add-ons/plugins for the acknowledgements section. Here are a few examples.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

This is an example of how you may give instructions on setting up your project locally.
To get a local copy up and running follow these simple example steps.

### Prerequisites

This is an example of how to list things you need to use the software and how to install them.
* npm
  ```sh
  npm install npm@latest -g
  ```

### Installation

_Below is an example of how you can instruct your audience on installing and setting up your app. This template doesn't rely on any external dependencies or services._

1. Get a free API Key at [https://example.com](https://example.com)
2. Clone the repo
   ```sh
   git clone https://github.com/github_username/repo_name.git
   ```
3. Install NPM packages
   ```sh
   npm install
   ```
4. Enter your API in `config.js`
   ```js
   const API_KEY = 'ENTER YOUR API';
   ```
5. Change git remote url to avoid accidental pushes to base project
   ```sh
   git remote set-url origin github_username/repo_name
   git remote -v # confirm the changes
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

- [x] Add Changelog
- [x] Add back to top links
- [ ] Add Additional Templates w/ Examples
- [ ] Add "components" document to easily copy & paste sections of the readme
- [ ] Multi-language Support
    - [ ] Chinese
    - [ ] Spanish

See the [open issues](https://github.com/othneildrew/Best-README-Template/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>
