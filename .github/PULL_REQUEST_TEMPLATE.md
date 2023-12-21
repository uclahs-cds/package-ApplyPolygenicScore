> This is a template for UCLA-CDS R package developers to create a github pull request template. Things should be adjusted for individual pipeline including:
> 1. additional checklist items specific to the package
> 2. a description of how testing is expected to be done
> 3. a template list or table for testing results
> 4. additional notes wrapped in \<!--- ---> (or \<!-- --\> for inline comments) that help PR submitters to fill in.
> 5. delete this block of instructional text.

<!--- Please read each of the following items and confirm by replacing the [ ] with a [ ] --->

- [ ] I have read the [code review guidelines](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3187646/Code+Review+Guidelines) and the [code review best practice on GitHub check-list](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3189956/Code+Review+Best+Practice+on+GitHub+-+Check+List).

- [ ] The name of the branch is meaningful and well formatted following the [standards](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3189956/Code+Review+Best+Practice+on+GitHub+-+Check+List), using \[AD_username (or 5 letters of AD if AD is too long)-\[brief_description_of_branch].

- [ ] I have set up or verified the branch protection rule following the [github standards](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3190380/GitHub+Standards#GitHubStandards-Branchprotectionrule) before opening this pull request.

- [ ] I have added the changes included in this pull request to `NEWS` under the next release version or unreleased, and updated the date.

- [ ] I have updated the version number in `metadata.yaml` and `DESCRIPTION`.

- [ ] Both `R CMD build` and `R CMD check` run successfully.

<!--- Briefly describe the changes included in this pull request and the test cases below
 !--- starting with 'Closes #...' if appropriate --->

Closes #...

## Testing Results

### Case 1
`input code`
`output`
### Case 2
`input code`
`output`