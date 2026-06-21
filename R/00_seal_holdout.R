# 00_seal_holdout.R  —  Day 1 of the preregistered confirmatory plan.
# Freezes, with a recorded random seed, an EXPLORE / CONFIRM split of each dataset.
# The preregistered confirmatory multiverse will run ONLY on the CONFIRM (held-out)
# rows. Run this ONCE; never edit the seed afterward.

SEED      <- 20260620L     # <-- record this exact number in the preregistration
CONF_FRAC <- 0.50          # fraction held out for confirmation (large N -> generous)

OUT <- "/Users/menghao/Library/Application Support/Claude/local-agent-mode-sessions/3aa1fdd3-f508-4003-af4a-2017e07db0cb/f260600f-e911-4881-9c9e-3e3f72bfd985/local_554904fd-182b-4ded-a87d-21199e51ca7f/outputs"
DEST <- "/Users/menghao/Downloads/study1_repro/prereg"

# read only the .npy header to get the row count (no array load)
npy_nrow <- function(path){
  con <- file(path, "rb"); on.exit(close(con))
  readBin(con, "raw", 8)                                   # magic + version
  hlen <- readBin(con, "integer", 1, size = 2, endian = "little")
  hdr  <- rawToChar(readBin(con, "raw", hlen))
  as.integer(gsub("[^0-9]", "", regmatches(hdr, regexpr("\\(([0-9]+),", hdr))))
}

sets <- list(NEO120 = "domain_scores.npy", IPIP50 = "ipip50_scores.npy", HEXACO = "hexaco_scores.npy")
explored_n <- 20000L       # exploratory pass used 20k subsamples (per the note)

set.seed(SEED)
lines <- c("HOLD-OUT MANIFEST  (Study 1: do personality types exist)",
           paste("sealed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
           paste("seed:", SEED, "| confirm fraction:", CONF_FRAC), "")
for (nm in names(sets)){
  N    <- npy_nrow(file.path(OUT, sets[[nm]]))
  perm <- sample.int(N)
  nc   <- floor(N * CONF_FRAC)
  confirm <- sort(perm[seq_len(nc)])                       # 1-based row indices
  writeLines(as.character(confirm), file.path(DEST, paste0("holdout_confirm_", nm, ".txt")))
  seen_frac <- round(100 * explored_n / N, 1)
  note <- if (seen_frac > 50) "  *** exploratory pass already touched most of this set; its hold-out is NOT naive -> treat HEXACO as robustness replication, not clean confirmation ***" else ""
  lines <- c(lines,
             sprintf("%-7s  total N = %7d  | explore = %7d  | CONFIRM(held-out) = %7d  (exploratory ~%s%% seen)%s",
                     nm, N, N - nc, nc, seen_frac, note))
}
lines <- c(lines, "",
  "USE: the confirmatory multiverse loads each .npy and subsets to the rows in",
  "holdout_confirm_<set>.txt (1-based; subtract 1 for Python). The EXPLORE rows",
  "(everything else) may be used for any further exploratory work.",
  "",
  "PREREG SENTENCE (paste): 'Prior to the confirmatory analysis we partitioned each",
  sprintf("dataset with a fixed seed (%d) into a 50%% exploration set and a 50%% held-out", SEED),
  "confirmation set (row indices frozen and timestamped on OSF); all confirmatory",
  "specifications were run only on the held-out set. Because the exploratory pass had",
  "already accessed ~20,000-row subsamples drawn before this split, we disclose this",
  "prior contact; for HEXACO, whose N is small, the hold-out cannot be treated as naive.'")
writeLines(lines, file.path(DEST, "HOLDOUT_MANIFEST.txt"))
cat(paste(lines, collapse = "\n"), "\n")
