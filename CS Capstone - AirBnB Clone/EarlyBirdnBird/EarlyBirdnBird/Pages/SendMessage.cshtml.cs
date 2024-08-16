using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using System.ComponentModel.DataAnnotations;

namespace EarlyBirdnBird.Pages
{
    public class SendMessageModel : PageModel
    {
        private readonly ApplicationDbContext _dbContext;
        private readonly UserManager<AppUser> _userManager;


        public SendMessageModel(ApplicationDbContext dbContext, UserManager<AppUser> userManager)
        {
            _dbContext = dbContext;
            _userManager = userManager;
        }

        public IActionResult OnGet()
        {
            return Page();
        }

        [BindProperty]
        public InputModel Input { get; set; }

        public class InputModel
        {
            [Required]
            [DataType(DataType.Text)]
            [Display(Name = "To")]
            public string To { get; set; }

            [Required]
            [DataType(DataType.Text)]
            [Display(Name = "Subject")]
            public string Subject { get; set; }

            [Required]
            [DataType(DataType.Text)]
            [Display(Name = "Message")]
            public string Body { get; set; }

        }

        public async Task<IActionResult> OnPostAsync()
        {

            //NEED ERROR HANDLING!!!!!!!!!!!!!!!!!!!

            var currentUser = await _userManager.GetUserAsync(User);

            if (currentUser == null)
            {
                // User not found, handle the error accordingly
                return NotFound("User not found.");
            }

            UserMail newMail = new UserMail()
            {
                fromUser = currentUser,
                toUser = Input.To,
                subject = Input.Subject,
                body = Input.Body,
                DateSent = DateTime.Now
            };

            _dbContext.UserMails.Add(newMail);
            await _dbContext.SaveChangesAsync();

            //return to inbox
            return RedirectToPage("./Inbox");
        }
    }
}
