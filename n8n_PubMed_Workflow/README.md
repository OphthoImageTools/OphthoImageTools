# PubMed Monthly Digest &mdash; n8n workflow

Importable [n8n](https://n8n.io) workflow that on the **25th of every month** runs a PubMed search across the major ophthalmology journals for **uveitis &amp; scleritis**, asks **Claude (Anthropic)** to summarize each abstract in 2&ndash;3 clinical sentences, and emails the digest as a styled HTML newsletter via Gmail.

## What the workflow does

```
Schedule (25th @ 9:00)
   │
   ▼
NCBI ESearch  ───────►  PubMed query (uveitis / scleritis ∩ ophthalmology journals)
   │
   ▼
Extract PMIDs (JS)
   │
   ▼
NCBI EFetch  ────────►  XML article details
   │
   ▼
Parse XML (JS)  ─────►  title, authors, journal, date, abstract, PMID
   │
   ▼
Claude Summarize  ───►  per-article 2&ndash;3 sentence clinical summary (JSON array)
   │
   ▼
Format Email (JS)  ──►  styled HTML body with article cards + PubMed links
   │
   ▼
Gmail Send  ─────────►  inbox
```

## File

| File | Role |
|---|---|
| `pubmed_workflow.json` | Importable n8n workflow JSON, sanitized of personal data and credentials. All sensitive values are replaced with `__REPLACE_ME_*__` placeholders that you fill in with your own values during setup. |

## Quick start (5 minutes)

### 1. Import the workflow

1. Download `pubmed_workflow.json` from this folder.
2. In n8n, click **Workflows** &rarr; ⋯ menu &rarr; **Import from File** &rarr; pick the JSON.
3. The workflow opens in a new tab named **PubMed Monthly Digest (OphthoImageTools template)**, set to **inactive** by default.

### 2. Replace the placeholders

Open each node and substitute the `__REPLACE_ME_*__` strings:

| Node | What to set |
|---|---|
| **Claude Summarize** | In the *Headers* parameter, replace `__REPLACE_ME_ANTHROPIC_API_KEY__` with your Anthropic API key (get one at https://console.anthropic.com/settings/keys). |
| **Send Gmail** | Replace `__REPLACE_ME_EMAIL_RECIPIENT__` with the email address you want the digest sent to. |
| **Send Gmail &rarr; Credentials** | Click *Create new credential* &rarr; **Gmail OAuth2** &rarr; complete the OAuth flow to authorize n8n to send mail from your account. |

> [!IMPORTANT]
> **Never paste an API key directly in a node parameter** as the original workflow did. The recommended pattern is to use n8n's **credential manager**: in the *Claude Summarize* node, switch *Authentication* to *Generic Credential Type* &rarr; *Header Auth* &rarr; create a new credential called e.g. "Anthropic API" with `x-api-key` as the name and your key as the value. Credentials live in n8n's encrypted store and are **never included in workflow exports**.

### 3. Customize the search (optional)

The default query in **PubMed Search** targets uveitis / scleritis across ~16 ophthalmology journals. You can edit it freely &mdash; the URL uses standard [PubMed E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) syntax. Examples:

- **Different topic**: replace the `uveitis[tiab] OR ...` block with your own boolean query.
- **Different journals**: add / remove `"Journal Name"[Journal]` entries.
- **Different time window**: the workflow uses *current month so far* (`mindate=YYYY/MM/01` &rarr; `maxdate=today`). Change to `relative date` (`reldate=30`) for "last 30 days regardless of month boundary".

Build and validate queries interactively at https://pubmed.ncbi.nlm.nih.gov/advanced/.

### 4. Customize the schedule (optional)

The trigger is a cron expression `0 9 25 * *` &mdash; **9:00 AM on the 25th of every month**. Edit it in the *Every 25th of Month* node:

- Weekly Monday: `0 9 * * 1`
- Daily 8 AM: `0 8 * * *`
- First Friday of the month: `0 9 * * 5#1` (depending on n8n version)

### 5. Activate

Click the **Active** toggle (top right of the workflow). The first run will fire at the next scheduled time. To test immediately, click **Execute Workflow** &mdash; you should receive the digest within ~30 seconds.

## Customizing the AI prompt

The clinical-summary tone is set inside the **Parse Articles** node (the `prompt` variable, at the bottom of the JS code). Default:

> *You are an expert ophthalmologist. For each article below, write a concise clinical summary in English (2&ndash;3 sentences) based on the abstract. Reply ONLY with a valid JSON array, no markdown, no backticks, no additional text. Format: [{"pmid":"...", "summary":"..."}]*

Tweak it for:
- **Different language** (Italian, Spanish, etc.).
- **Different audience** (e.g. *summarize for residents* vs. *for senior clinicians*).
- **Highlight specific elements** (study design, sample size, key biomarkers, treatment effect size).

Keep the *"Reply ONLY with a valid JSON array..."* part &mdash; the next node (`Format Email`) parses the response as JSON.

## Cost notes

- **PubMed E-utilities**: free, but rate-limited. With 100+ articles per month you should register an [NCBI API key](https://account.ncbi.nlm.nih.gov) and add `&api_key=YOUR_KEY` to both NCBI URLs.
- **Anthropic Claude**: token cost depends on volume. ~100 abstracts &rarr; ~30k input tokens + ~10k output tokens. With Claude Sonnet at current pricing this is roughly **$0.10 per monthly run**.
- **Gmail**: free; Gmail's daily send limit (500 / 2,000 depending on account type) is irrelevant for a personal digest.

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| Workflow runs but no email arrives | Gmail OAuth credential expired &rarr; re-authorize in *Credentials*. |
| `401 Unauthorized` on Claude Summarize | API key wrong or revoked. Generate a new one and update the credential. |
| Digest body shows "Summary not available" for every article | Claude returned malformed JSON. Open the node, switch to *Execute Once* mode, inspect the raw response. Usually a model change or a prompt drift &mdash; tighten the *"Reply ONLY with a valid JSON array..."* instruction. |
| ESearch returns 0 articles | Query syntax error or month boundary issue. Test the URL by pasting it in a browser. |

## License &amp; attribution

This workflow is part of [OphthoImageTools](https://github.com/OphthoImageTools/OphthoImageTools), MIT-licensed. If you adapt it for a publication or talk, a citation to the parent repository is appreciated &mdash; see the root [`CITATION.cff`](../CITATION.cff).
