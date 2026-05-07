# PubMed automation &mdash; n8n workflow

Importable [n8n](https://n8n.io) workflow that runs scheduled PubMed searches on topics of interest and emails / chats a curated digest.

## Status

The cleaned, sharable workflow JSON will be added to this folder shortly. The export will be **stripped of personal credentials** (NCBI API key, email recipient, OpenAI key etc.) and replaced with placeholders that you can fill in with your own values.

## What the workflow does (overview)

1. **Schedule trigger** &mdash; runs daily / weekly.
2. **NCBI E-utilities &mdash; ESearch** &mdash; queries PubMed for the user-defined topic (e.g. `"uveitis"[Title/Abstract] AND ("OCT angiography" OR "macrophage-like cells")`).
3. **NCBI E-utilities &mdash; EFetch** &mdash; pulls the abstract, authors, journal, DOI for each PMID.
4. **De-duplication** against the previously-stored PMID list.
5. **(Optional)** LLM summarization step that produces a short Italian / English bullet-summary per article.
6. **Output** &mdash; email or messaging-app digest, plus persistence to a sheet / Notion / database.

## How to import (once the JSON is uploaded)

1. Download `pubmed_workflow.json` from this folder.
2. In n8n, click **Workflows &rarr; Import from File**.
3. Open the imported workflow and replace every node where credentials are marked `__REPLACE_ME__`:
   - **NCBI API key** &mdash; get one free at https://account.ncbi.nlm.nih.gov &rarr; *API Key Management*.
   - **Email / Slack / Telegram credentials** &mdash; whichever output you choose.
   - **(Optional)** OpenAI/Anthropic key for the summary node.
4. Edit the *Schedule* trigger to your preferred cadence.
5. Edit the *Search query* in the ESearch node to match your interests (PubMed advanced search syntax is supported).
6. Click **Save** &rarr; **Activate**.

## Customizing the search

The default query is a placeholder. Use [PubMed advanced search](https://pubmed.ncbi.nlm.nih.gov/advanced/) to build your own and paste it into the `term` parameter of the ESearch node. Examples:

- `"vogt-koyanagi-harada"[Title/Abstract] AND 2024:2026[dp]`
- `("keratic precipitates"[tiab] OR "anterior chamber inflammation"[tiab]) AND ("optical coherence tomography"[mh] OR OCT[tiab])`
- `(uveitis[mh] OR uveitis[tiab]) AND (artificial intelligence[tiab] OR machine learning[tiab] OR deep learning[tiab])`

## To do

- [ ] Upload the cleaned `pubmed_workflow.json`.
- [ ] Add screenshots of the workflow graph.
- [ ] Add a step-by-step setup video link.
