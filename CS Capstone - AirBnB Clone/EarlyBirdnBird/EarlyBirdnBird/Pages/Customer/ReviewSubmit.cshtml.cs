using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Microsoft.EntityFrameworkCore;
using System.ComponentModel.DataAnnotations;

namespace EarlyBirdnBird.Pages.Customer
{
    public class ReviewSubmitModel : PageModel
    {
        private readonly ApplicationDbContext _dbContext;
        private readonly UserManager<AppUser> _userManager;

        public ReviewSubmitModel(ApplicationDbContext dbContext, UserManager<AppUser> userManager)
        {
            _dbContext = dbContext;
            _userManager = userManager;
        }

        public Review Reviews { get; set; }
        public Rental Rentals { get; set; }
        private int reviewId;

        [BindProperty]
        public InputModel Input { get; set; }

        public class InputModel
        {
            [Required]
            [DataType(DataType.Text)]
            [Display(Name = "Message")]
            public string Body { get; set; }
        }

        public async Task OnGetAsync(int? id)
        {
            Rentals = await _dbContext.Rentals.FirstOrDefaultAsync(r => r.Id == id);
        }

        public async Task<IActionResult> OnPostAsync(int? id)
        {

            //NEED ERROR HANDLING!!!!!!!!!!!!!!!!!!!

            var currentUser = await _userManager.GetUserAsync(User);

            reviewId = (int)id;

            if (currentUser == null)
            {
                // User not found, handle the error accordingly
                return NotFound("User not found.");
            }

            Review newReview = new Review()
            {
                RentalId = reviewId,
                Reviewer = currentUser,
                body = Input.Body
            };

            _dbContext.Reviews.Add(newReview);
            await _dbContext.SaveChangesAsync();

            //return to home
            return RedirectToPage("./Home");
        }
    }
}
